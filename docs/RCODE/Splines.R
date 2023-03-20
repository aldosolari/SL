#=======================================
# SPLINES
#=======================================

#---------------------------------------
# Example
#---------------------------------------

rm(list=ls())

n <- 500
sigma <- 0.3
set.seed(123)
x = sort(runif(n)) 
f_x = sin(2*(4*x-2)) + 2*exp(-(16^2)*((x-.5)^2))
y = f_x + rnorm(n, mean=0, sd=sigma)

#---------------------------------------
# Smoothing spline
#---------------------------------------

#pdf("Figure_smooth_spline.pdf")
overfit <- smooth.spline(x, y, all.knots=T, spar = 0)
fit_smooth <- smooth.spline(x, y, all.knots=T, cv=TRUE)
fit_smooth
plot(x,y,col="gray")
lines(x, overfit$y, col="gray")
lines(x, f_x, col=2, lwd=2)
lines(x, fit_smooth$y, col=4, lwd=2)
#dev.off()

#---------------------------------------
# B-splines
#---------------------------------------

tpower <- function(x, t, deg){
  (x - t) ^ deg * (x > t)
}

bbase <- function(x, xl, xr, ndx, deg){
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  Delta <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(Delta)
  B
}


xl=min(x)
xr=max(x)
ndx=10
bdeg=3

B <- bbase(x, xl, xr, ndx, bdeg)
knots_all <- seq(xl - bdeg * (xr - xl) / ndx, xr + bdeg * (xr - xl) / ndx, by = (xr - xl) / ndx)
#pdf("Figure_bspline.pdf")
plot(knots_all,rep(0,length(knots_all)),pch=19, ylab=expression(B[j](x)), xlab="x")
abline(v=knots_all, lty=2)
for (i in 1:ncol(B)) lines(x,B[,i], col=i, lwd=2)
#dev.off()

rowSums(B)

max(svd(B)$d) / min(svd(B)$d)

y_hat <- B %*% solve(crossprod(B)) %*% crossprod(B, y)

#---------------------------------------
# P-splines
#---------------------------------------

lambda <- 0.18
O <- 2

D <- diag(ncol(B))
for (k in 1:O) D <- diff(D)
beta_hat <- solve(t(B) %*% B + lambda * t(D) %*% D, t(B) %*% y)
y_hat <- B %*% beta_hat

RSS <- sum((y - y_hat)^2)
tr_S <- sum(diag(solve(t(B) %*% B + lambda * t(D) %*% D) %*% (t(B) %*% B)))
GCV <- (1/n) * ( RSS / ( 1 - tr_S / nrow(B) )^2 )


plot(x,y ,col="gray")
lines(x, f_x, col=2, lwd=2)
abline(v=knots_all, lty=2)
lines(x,y_hat, lwd=2)

library(JOPS)
fit <- psNormal(x, y, nseg = 10, bdeg = 3, pord = 2, lambda=lambda)
sum(abs(fit$muhat - y_hat))

#---------------------------------------
# mcycle data 
#---------------------------------------

rm(list=ls())

data(mcycle)
x = mcycle$times
y = mcycle$accel

