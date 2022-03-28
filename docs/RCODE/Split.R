#=======================================
# DATA SPLITTING FOR VARIABLE SELECTION
#=======================================

#---------------------------------------
# NAIVE TWO-STEP PROCEDURE
#---------------------------------------

rm(list=ls())

alpha <- 0.05
n <- 200
p <- 1000
s <- 10   

beta <- c(rep(1,s),rep(0,p-s))
S <- which(beta != 0)
varType <- rep("N",p)			
varType[S] <- "S"
rho <- 0
Sigma <- toeplitz(rho^(0:(p-1)))
SNR <- 2.5
sigma2 <- (t(beta) %*% Sigma %*% beta) / SNR

set.seed(123)
X <- as.matrix(matrix(rnorm(n * p), n, p) %*% chol(Sigma) )
y <- X %*% beta + rnorm(n,mean=0,sd=sqrt(sigma2))

fit <- cv.glmnet(X, y)
M_hat <- which(coef(fit, s=fit$lambda.1se)[-1] != 0)
table(varType[M_hat])
m_hat <- length(M_hat)
M_hat_typeI <- length(setdiff(M_hat,S))
M_hat_typeII <- length(setdiff(S,M_hat))

fit_M_hat <- lm(y ~ X[,M_hat])
pval_M_hat <- summary(fit_M_hat)$coef[-1,4]	
S_hat <- M_hat[(pval_M_hat <= alpha)]
table(varType[S_hat])
s_hat = length(S_hat)
S_hat_typeI <- length(setdiff(S_hat,S))
S_hat_typeII <- length(setdiff(S,S_hat))

#---------------------------------------
# SINGLE-SPLIT
#---------------------------------------

set.seed(123)
L <- as.logical(sample(rep(0:1, each=n/2)))
I = !L
fit_L <- cv.glmnet(X[L,], y[L])
M_hat <- which(coef(fit_L, s=fit_L$lambda.1se)[-1]!=0)
table(varType[M_hat])

fit_I <- lm(y[I]~X[I, M_hat])
pval = rep(1,p)
pval[M_hat] = summary(fit_I)$coefficients[-1,4]
S_hat <- M_hat[pval[M_hat] <= alpha]
table(varType[S_hat])

pval_tilde <- rep(1,p)
pval_tilde[M_hat] <- p.adjust(pval[M_hat],"bonf")
S_tilde <- M_hat[pval_tilde[M_hat] <= alpha]
table(varType[S_tilde])

#---------------------------------------
# P-VALUE LOTTERY
#---------------------------------------

B <- 25
pval_matrix <- matrix(1,ncol=p,nrow=B)
pval_matrix_tilde <- pval_matrix

set.seed(123)
for (i in 1:B) {
  split <- as.logical(sample(rep(0:1, each=n/2)))
  fit <- cv.glmnet(X[split,], y[split])
  M_hat <- which( coef(fit, s=fit$lambda.1se)[-1] != 0 )
  fit <- lm(y[!split]~X[!split, M_hat])
  pval_matrix[i, M_hat] <- summary(fit)$coeff[-1,4]
  pval_matrix_tilde[i, M_hat] <- p.adjust(pval_matrix[i, M_hat], "holm") 
}
#pdf("Figure_plottery.pdf")
hist(pval_matrix[,S[1]], main=paste(B,"random splits"), xlab="p-value", 20)
#dev.off()

#pdf("Figure_pmedian.pdf")
boxplot(pval_matrix_tilde[,S], label=S)
abline(h=alpha/2, lty=3)
#dev.off()
pval_aggr <- pmin(2*apply(pval_matrix_tilde,2,median),1)
sum(pval_aggr <= alpha)

#---------------------------------------
# MULTI-SPLIT
#---------------------------------------

B = 25
library(hdi)
set.seed(123)
fit <- multi.split(x=X, y=y, B=B, fraction=0.5,
                   model.selector = lasso.cv,	 
                   ci = TRUE, ci.level = 1- alpha, 
                   gamma = 0.5)
S_hat <- which(fit$pval.corr <= alpha)
table(varType[S_hat])
confint(fit)[S_hat,]

#---------------------------------------
# SIMULATION
#---------------------------------------

sim_naive <- function(n,p,s,SNR,rho=0,alpha=0.05){
  Sigma <- toeplitz(rho^(0:(p-1)))
  beta = c(rep(1,s),rep(0,p-s))
  S <- which(beta != 0)
  sigma2 <- (t(beta) %*% Sigma %*% beta) / SNR
  X <- as.matrix(matrix(rnorm(n * p), n, p) %*% chol(Sigma) )
  y <- X %*% beta + rnorm(n,mean=0,sd=sqrt(sigma2))
  fit <- cv.glmnet(X, y)
  M_hat <- which(coef(fit, s=fit$lambda.1se)[-1] != 0)
  s_hat <- 0
  S_hat_typeI = 0
  S_hat_typeII  = s
  if( length(M_hat) > 0){
      fit_M_hat <- lm(y ~ X[,M_hat])
      pval_M_hat <- summary(fit_M_hat)$coef[-1,4]	
      S_hat <- M_hat[(pval_M_hat <= alpha)]
      s_hat = length(S_hat)
      S_hat_typeI <- length(setdiff(S_hat,S))
      S_hat_typeII <- length(setdiff(S,S_hat))
  }
  return(c(s_hat = s_hat, typeI=S_hat_typeI, typeII = S_hat_typeII))
}

B <- 10
my_n <- 200
my_p <- 1000
my_s <- 10
my_SNR <- 2.5
my_rho = 0
my_alpha <- 0.05
set.seed(123)
res = replicate(B, sim_naive(n = my_n, p = my_p, s = my_s, SNR = my_SNR, rho = my_rho, alpha=my_alpha))
mean(res[2,]>0)
mean(apply(res,2, function(x) ifelse(x[1]>0, x[2]/x[1],0)))
mean( (res[1,]-res[2,])/my_s )

