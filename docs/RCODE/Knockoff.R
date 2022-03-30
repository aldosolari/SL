#=======================================
# KNOCKOFF FILTER
#=======================================

rm(list=ls())

n <- 1000
p <- 200
s <- 40   
set.seed(123)
beta <- sample(c(rep(2,s/2),rep(-2,s/2),rep(0,p-s)))
S <- which(beta != 0)
N <- setdiff(1:p,S)
varType <- rep("N",p)			
varType[S] <- "S"
rho <- 0
Sigma <- toeplitz(rho^(0:(p-1)))
sigma2 <- 1

normc = function(X,center=T) {
  X.centered = scale(X, center=center, scale=F)
  X.scaled = scale(X.centered, center=F, scale=sqrt(colSums(X.centered^2)))
  X.scaled[,]
}

X_raw <- as.matrix(matrix(rnorm(n * p), n, p) %*% chol(Sigma) )
X <- normc(X_raw, center=T)
y_raw <- X %*% beta + rnorm(n,mean=0,sd=sqrt(sigma2))
y <- normc(y_raw, center=T)

#---------------------------------------
# KNOCKOFF CONSTRUCTION
#---------------------------------------

X_svd = svd(X)
Q = qr.Q(qr(cbind(X_svd$u, matrix(0,n,p))))
U = Q[,(p+1):(2*p)]

Sigma_inv = solve(crossprod(X))
d = 0.6
CtC = 4*( (d/2) * diag(rep(1,p)) - (d/2)^2 * Sigma_inv)
C = chol(CtC)
Xtilde = X %*% ( diag(rep(1,p)) - d * Sigma_inv  ) + U %*% C

crossprod(X)[1:3,1:3]
crossprod(Xtilde)[1:3,1:3]
crossprod(X,Xtilde)[1:3,1:3]

#library(knockoff)
#Xtilde = create.fixed(X,method="sdp")$Xk

#---------------------------------------
# KNOCKOFF STATISTICS
#---------------------------------------

XX <- cbind(X,Xtilde)
fit <- glmnet(XX, y)
l <- cv.glmnet(XX, y)$lambda.min
varType <- c(rep("N",p),rep("K",p)) 			
varType[S] <- "S"

S_hat <- which(coef(fit, s=l)[-1] != 0)
table(varType[S_hat])

#pdf("Figure_lasso_knockoff_true_FDP.pdf")
#op <- par(mar=c(4,5,4,4))
plot(1:(2*p),coef(fit, s=l)[-1], xlab="Index j", ylab=expression(hat(beta)[j]), pch=19, col=1+2*(varType=="K")+(varType!="N"))
#par(op)
#dev.off()

first_nonzero <- function(x) match(T, abs(x) > 0)
indices <- apply(fit$beta, 1, first_nonzero)
Z = ifelse(is.na(indices), 0, fit$lambda[indices] * n)
orig = 1:p
W = pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])

#---------------------------------------
# KNOCKOFF FDP ESTIMATE
#---------------------------------------

tau = 2

#pdf("Figure_knockoff_stat.pdf")
#op <- par(mar=c(4,5,4,4))
plot(Z[orig], Z[orig+p], col=1+(varType=="S"), pch=19, asp=1, 
     xlab=expression(paste(lambda," when ",X[j]," enters ")),
     ylab=expression(paste(lambda," when ",tilde(X)[j]," enters ")))
abline(a=0,b=1)
abline(h=tau,lty=3)
abline(v=tau,lty=3)
#par(op)
#dev.off()

S_hat_tau = which(W >= tau)
table(varType[S_hat_tau])
(1 + sum(W <= -tau)) / length(S_hat_tau)
sum(W[N] >= tau) / sum(W >= tau)


taus = sort(c(0,abs(W)))

FDP_hat = sapply(taus, function(tau)
  (1 + sum(W <= -tau)) / max(1, sum(W >= tau)))

FDP_true = sapply(taus, function(tau)
  (sum(W[varType=="N"] >= tau)) / max(1, sum(W >= tau)))

#pdf("Figure_knockoff_FDP.pdf")
plot(taus,FDP_true, type="l", lwd=2, xlab=expression(tau), ylab="FDP")
lines(taus,FDP_hat, col=2, lwd=2)
legend("topright", col=1:2, c("True","Estimate"), lty=1)
#dev.off()

alpha = 0.1
tau_hat <- taus[which(FDP_hat  <= alpha)[1]]
S_hat = which(W >= tau_hat)
table(varType[S_hat])
(1 + sum(W <= -tau_hat)) / length(S_hat)
sum(W[N] >= tau_hat) / sum(W >= tau_hat)


#---------------------------------------
# VARIABLE IMPORTANCE STATISTICS
#---------------------------------------

library(ranger)

random_forest_importance <- function(X, y, ...) {
  df = data.frame(y=y, X=X)
  rfFit = ranger::ranger(y~., data=df, importance="impurity", write.forest=F, ...)
  as.vector(rfFit$variable.importance)
}

set.seed(123)

Z = random_forest_importance(cbind(X, Xtilde), y) 
W = abs(Z[orig]) - abs(Z[orig+p])

tau = 0.001

#pdf("Figure_varimp.pdf")
plot(W, col=1+(varType=="S"), type="h", lwd=2, xlab="Index j")
abline(h=tau,lty=3)
abline(h=-tau,lty=3)
#dev.off()

S_hat_tau = which(W >= tau)
table(varType[S_hat_tau])
(1 + sum(W <= -tau)) / length(S_hat_tau)
sum(W[N] >= tau) / sum(W >= tau)


taus = sort(c(0, abs(W)))

FDP_hat = sapply(taus, function(tau)
  (1 + sum(W <= -tau)) / max(1, sum(W >= tau)))

FDP_true = sapply(taus, function(tau)
  (sum(W[N] >= tau)) / max(1, sum(W >= tau)))

#pdf("Figure_varimp_FDP.pdf")
plot(taus,FDP_true, type="l", lwd=2, xlab=expression(tau), ylab="FDP", ylim=c(0,1))
lines(taus,FDP_hat, col=2, lwd=2)
legend("topright", col=1:2, c("True","Estimate"), lty=1)
#dev.off()

alpha = 0.1
tau_hat <- taus[which(FDP_hat  <= alpha)[1]]
S_hat = which(W >= tau_hat)
table(varType[S_hat])
(1 + sum(W <= -tau_hat)) / length(S_hat)
sum(W[N] >= tau_hat) / sum(W >= tau_hat)

#---------------------------------------
# MODEL-X KNOCKOFF
#---------------------------------------

library(knockoff)

set.seed(123)

mu = rep(0,p)
gaussian_knockoffs = function(X) create.gaussian(X, mu, Sigma)
result = knockoff.filter(X_raw, y_raw, knockoffs=gaussian_knockoffs)
print(result)

fdp = function(selected) sum(beta[selected] == 0) / max(1, length(selected))
fdp(result$selected)
