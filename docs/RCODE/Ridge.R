#=======================================
# RIDGE REGRESSION
#=======================================

#---------------------------------------
# TOY EXAMPLE
#---------------------------------------

rm(list=ls())

X <- matrix(c(10^9, -1, -1, 10^(-5)), 2, 2)
beta <- c(1,1)
y <- X %*% beta

d <- svd(crossprod(X))$d
max(d)/min(d)
kappa(crossprod(X), exact=T)
solve( crossprod(X), crossprod(X, y) )
.Machine$double.eps

1/kappa(X, exact=T)
solve(X,y)

#---------------------------------------
# CEMENT DATA
#---------------------------------------

rm(list=ls())

library(MASS)
library(glmnet)
library(genridge)

y <- cement$y
X <- data.matrix(cement[,-5])
p <- ncol(X)
n <- nrow(X)

cor(X)
fit <- lm(y ~., cement)
summary(fit)
car::vif(fit)

l = c(glmnet(X,y,alpha=0)$lambda, seq(1,0,length.out = 100))
fit_ridge <- lm.ridge(y ~., cement, lambda = l)

#pdf("Figure_ridge_trace.pdf")
plot(log1p(l),coef(fit_ridge)[,2], 
     xlim=range(log1p(l)), ylim=c(-0.5,1.5), type="l", ylab="Coefficient", xlab=expression(log(lambda+1)))
for (i in 1:p) lines(log1p(l), coef(fit_ridge)[,i+1], type="l", col=i)
points(rep(0,4),fit$coefficients[-1], col=1:p)
abline(h=0, lty=3)
#dev.off()


#pdf("Figure_ridge_pairs.pdf")
pairs(ridge(y, X, lambda=c(0, 0.1, 1, 10, 1000)), radius=0.5)
#dev.off()

#---------------------------------------
# OVERFITTING 
#---------------------------------------

rm(list=ls())

y = c(-1,0)
X = matrix(c(.75,.5, 1,.5),ncol=2)
coef(lm(y~0+X))

#pdf("Figure_overfitting.pdf")
plot(y[1],y[2], xlim=c(-3,3), ylim=c(-1,3), asp=1, xlab=expression(u[1]), ylab=expression(u[2]))
arrows(x0=0,y0=0, x1=y[1], y1=y[2], length = 0.1, lwd=2)
text(y[1],y[2], expression(y), pos=1)
text(2.5,1,expression(paste(beta[1], "= 4")), col=2)
text(2.5,0.5,expression(paste(beta[2], "= -4")), col=4)
abline(h=0)
abline(v=0)
arrows(x0=0,y0=0, x1=X[1,1], y1=X[2,1], length = 0.1, col=2, lwd=3)
text(X[1,1], X[2,1], expression(x[1]), pos=3)
for (k in 2:4) arrows(x0=0,y0=0, x1=k*X[1,1], y1=k*X[2,1], length = 0.1 , col=2)
arrows(x0=0,y0=0, x1=X[1,2], y1=X[2,2], length = 0.1, col=4, lwd=3)
text(X[1,2], X[2,2], expression(x[2]), pos=4)
for (k in 1:4) arrows(x0=4*X[1,1],y0=4*X[2,1], x1=4*X[1,1]-k*X[1,2], y1=4*X[2,1]-k*X[2,2], length = 0.1 , col=4)
#dev.off()

#---------------------------------------
# MY RIDGE 
#---------------------------------------

rm(list=ls())

library(MASS)

my_ridge <- function(X, y, lambda){
  n <- nrow(X)
  p <- ncol(X)
  y_mean <- mean(y)
  y <- y - y_mean
  X_mean <- colMeans(X)
  X <- X - rep(1,n) %*% t(X_mean)
  X_scale <- sqrt( diag( (1/n) * crossprod(X) ) )
  X <- X %*% diag( 1 / X_scale )
  beta_scaled <- solve(crossprod(X) + lambda*diag(rep(1,p)), t(X) %*% y) 
  beta <- diag( 1 / X_scale ) %*% beta_scaled
  beta0 <- y_mean - X_mean %*% beta
  return(c(beta0, beta))
}

y <- cement$y
X <- data.matrix(cement[,-5])
n <- nrow(X)

l = 1
my_ridge(X,y,lambda = l)
lm.ridge(y ~ ., cement, lambda = l)
coef(glmnet(X, y, alpha=0, lambda = l/n, thresh = 1e-20))

y_std <- scale(y, center=TRUE, scale=sd(y)*sqrt((n-1)/n) )[,]
my_ridge(X,y_std,lambda = l)
coef(glmnet(X, y_std, alpha=0, lambda = l/n, thresh = 1e-20))


#---------------------------------------
# MSE
#---------------------------------------

rm(list=ls())

library(readr)

dataset <- read_csv("https://hastie.su.domains/CASI_files/DATA/diabetes.csv")[,-1]
X <- data.matrix(dataset[,-11])
y <- dataset$prog
n <- nrow(X)
p <- ncol(X)
Z <- scale(X, center=T, scale=sqrt(diag(var(X)*(n-1)/n)))[,]
beta <- matrix(coef(lm(I(y-mean(y)) ~ 0+Z)),ncol=1)
sigma2 <- summary(lm(I(y-mean(y)) ~ 0+Z))$sigma^2

ridge_MSE <- function(X,beta,sigma2,lambda){
  n <- nrow(X)
  p <- ncol(X)
  beta <- matrix(beta,ncol=1)
  SVD <- svd(X)
  d <- SVD$d
  U <- SVD$u
  V <- SVD$v
  Bias <- V %*%  diag( lambda/(d^2+lambda) ) %*% t(V) %*% beta 
  Var <- sigma2 * V %*% diag( ((d^2)/(d^2+lambda)^2) ) %*% t(V) 
  MSE <- sum(diag(Var)) + crossprod(Bias)
  return(c(MSE, crossprod(Bias), sum(diag(Var)) ) )
}

lambdas <- seq(0,5,length.out=100)
MSE <- sapply(lambdas, function(l) ridge_MSE(Z, beta, sigma2, lambda=l))

#pdf("Figure_ridge_MSE.pdf")
plot(lambdas, MSE[1,], xlab=expression(lambda), ylab="MSE", type="l", ylim=c(0,max(MSE[1,])), lwd=2)
abline(h=MSE[1,1], lty=3)
lines(lambdas, MSE[2,], col=2)
lines(lambdas, MSE[3,], col=3)
legend("bottomright", c("MSE","Bias2","Var"), col=1:3, lty=1)
#dev.off()


#---------------------------------------
# EPE
#---------------------------------------

rm(list=ls())

library(glmnet)

y <- longley[, "Employed"]
X <- data.matrix(longley[, c(2:6,1)])
n <- nrow(X)
p <- ncol(X)
Z <- scale(X, center=T, scale=sqrt(diag(var(X)*(n-1)/n)))[,]
beta <- matrix(coef(lm(I(y-mean(y)) ~ 0+Z)),ncol=1)
sigma2 <- summary(lm(I(y-mean(y)) ~ 0+Z))$sigma^2

ridge_EPE <- function(X,beta,sigma2,lambda){
  n <- nrow(X)
  p <- ncol(X)
  beta <- matrix(beta,ncol=1)
  SVD <- svd(X)
  d <- SVD$d
  U <- SVD$u
  V <- SVD$v
  Bias <- V %*%  diag( lambda/(d^2+lambda) ) %*% t(V) %*% beta 
  Var <- sigma2 * V %*% diag( ((d^2)/(d^2+lambda)^2) ) %*% t(V) 
  EPE <- mean(apply(X,1,function(x) t(x)%*%Var%*%x + (x%*% Bias)^2 )) + sigma2
  return(EPE)
}

l <- seq(0,.05,length.out=100)

set.seed(123)
y <- Z %*% beta + rnorm(n, 0, sd=sqrt(sigma2))
y <- y - mean(y)
cv_fit <- cv.glmnet(Z, y, alpha = 0, standardize = F, nfolds=n, grouped=FALSE, lambda = l)
l <- cv_fit$lambda
fit <- glmnet(Z, y, alpha = 0, standardize = F, lambda = l)

EPE <- sapply(l, function(l) ridge_EPE(Z, beta, sigma2, lambda=l))
LOO <- cv_fit$cvm

#pdf("Figure_ridge_EPE.pdf")
plot(l, EPE, xlab=expression(lambda), ylab="Prediction error", type="l",  lwd=2, ylim=c(min(EPE), max(LOO)))
lines(l,LOO)
legend("bottomright", c("LOO", "EPE"), lwd=1:2)
abline(v=l[which.min(EPE)], lty=3)
abline(v=l[which.min(LOO)], lty=2)
points(l[which.min(LOO)], LOO[which.min(LOO)])
points(l[which.min(EPE)], EPE[which.min(EPE)], pch=19)
#dev.off()

#pdf("Figure_ridge_Longley.pdf")
matplot(l, t(coef(fit)[-1,]), type = "l", lty=1, 
        xlab=expression(lambda), ylab=expression(hat(beta)[lambda]), 
        col=1:p, ylim=range(beta)) 
abline(v=l[which.min(EPE)], lty=3)
points(rep(0,p), coef(fit)[-1,length(l)], col=1:p)
points(rep(-.0015,p), beta, col=1:p, pch=19)
#dev.off()

#---------------------------------------
# CROSS-VALIDATION
#---------------------------------------

rm(list=ls())

library(readr)
library(glmnet)

dataset <- read_csv("https://hastie.su.domains/CASI_files/DATA/diabetes.csv")[,-1]
X <- data.matrix(dataset[,-11])
y <- dataset$prog
n <- nrow(X)
p <- ncol(X)
Z <- scale(X, center=T, scale=sqrt(n*diag(var(X)*(n-1)/n)))[,]
colnames(Z)<- names(dataset)[-1]

cv_fit <- cv.glmnet(Z, y, alpha = 0, standardize = FALSE, nfolds=n, grouped=FALSE)
#pdf("Figure_ridge_cv.pdf")
plot(cv_fit)
#dev.off()
l_min <- cv_fit$lambda.min

l <- seq(0,0.25, length.out = 100)
fit <- glmnet(Z, y, alpha = 0, family="gaussian", standardize = FALSE, lambda = l)

#pdf("Figure_ridge_Diabetes.pdf")
matplot(l, t(coef(fit)[-1,length(l):1]), type = "l", lty=1, 
        xlab=expression(lambda), ylab=expression(hat(beta)[lambda]), 
        xlim = c(-.01, 0.25), col=1:p) 
text(x=-.01, y=coef(fit)[-1,length(l)], labels =names(dataset)[-1], cex=0.5)
abline(v=l_min, lty=3)
#dev.off()

#---------------------------------------
# KERNEL TRICK
#---------------------------------------

rm(list=ls())

library(readr)
dataset <- read_csv("https://hastie.su.domains/CASI_files/DATA/diabetes.csv")[,-1]
X_raw <- data.matrix(dataset[,-11])
y <- dataset$prog
n <- nrow(X_raw)
p <- ncol(X_raw)
X <- scale(X_raw, center=T, scale=sqrt(diag(var(X_raw)*(n-1)/n)))[,]

lambda = 1

X2 = cbind(X, do.call(cbind, lapply(1:p, function(i) X[,i] *  apply(X,2,identity)) ))
K2 <- X2 %*% t(X2)
yhat2 = K2 %*% solve(K2 + lambda*diag(n)) %*% y

K = sapply(1:n, function(j)
  apply(X,1,function(x) (0.5 + t(x) %*% X[j,] )^2 - 0.25 )
)
yhat = K %*% solve(K + lambda*diag(n)) %*% y




