#=======================================
# STABILITY SELECTION
#=======================================

rm(list=ls())

n <- 200
p <- 1000
s <- 10   

beta <- c(rep(1,s),rep(0,p-s))
S <- which(beta != 0)
varType <- rep("N",p)			
varType[S] <- "S"
rho <- 0
Sigma <- toeplitz(rho^(0:(p-1)))
SNR <- 1
sigma2 <- (t(beta) %*% Sigma %*% beta) / SNR

set.seed(123)
X <- as.matrix(matrix(rnorm(n * p), n, p) %*% chol(Sigma) )
y <- X %*% beta + rnorm(n,mean=0,sd=sqrt(sigma2))

# REGULARIZATION PATH--------------------------

fit <- glmnet(X, y)
Lambda <- fit$lambda
l <- cv.glmnet(X, y)$lambda.1se
S_hat <- which(coef(fit, s=l)[-1] != 0)
table(varType[S_hat])


col <- rep("lightgray", p)
col[S]<-"red"
#pdf("Figure_reg_path.pdf")
plot(fit, xvar="lambda", col=col, lwd=2)
#dev.off()
abline(v=log(l))

# STABILITY PATH--------------------------------

B = 100
S_hat_half <- array(NA, dim=c(B, p, length(Lambda)), dimnames=list(1:B, colnames(X), Lambda))

for (b in 1:B) {
  ind <- as.logical(sample(rep(0:1, each=n/2)))
  fit_b <- glmnet(X[ind,], y[ind], lambda=Lambda)
  S_hat_half[b,,] <- as.matrix(coef(fit_b)[-1,]!=0)
}
pi_hat <- apply(S_hat_half, 2:3, mean)
#pdf("Figure_stab_path.pdf")
matplot(log(Lambda), t(pi_hat), 
        type="l", lty=1, xlim=log(range(Lambda)), col=col, lwd=2, las=1, bty="n", xlab="Log Lambda", ylab="Estimated probability", ylim=c(0,1))
#dev.off()
tau = 0.6
abline(h=tau)

library(stabs)
fit <- stabsel(x = X, y = y, fitfun = glmnet.lasso, q = 50,
               cutoff=0.6, assumption ="none") 
#pdf("Figure_stab.pdf")
plot(1:p,fit$max, xlim=c(1,p),col=col, xlab="Variable j", ylab=expression(hat(pi)[j]), pch=19)
abline(h=0.6, lty=3)
#dev.off()
plot(fit)
