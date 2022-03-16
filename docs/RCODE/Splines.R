#=======================================
# SPLINES
#=======================================


#---------------------------------------
# EXAMPLE
#---------------------------------------

rm(list=ls())

n <- 500
sigma <- 0.3
set.seed(123)
x = sort(runif(n)) 
f_x = sin(2*(4*x-2)) + 2*exp(-(16^2)*((x-.5)^2))
y = f_x + rnorm(n, mean=0, sd=sigma)

#pdf("Figure_nonlinear.pdf")
plot(x, y, col="gray") 
lines(x, f_x, col=2, lwd=2)
lines(x, fitted(lm(y ~ poly(x,15))), lwd=2)
legend("bottomright", col=2:1, c("f(x)",expression(hat(f)(x))), lty=1)
#dev.off()


#pdf("Figure_hatmatrix.pdf")
X = model.matrix(lm(y ~ poly(x,degree=15)))
S = X %*% solve(crossprod(X)) %*% t(X)
filled.contour(apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.001,0.001,max(S)), col=c("red", "white", "blue"))
#dev.off()

M = 3
knots = seq(0.1, 0.9, by=.1)
x_cut = cut(x, c(-Inf,knots,Inf) )
#pdf("Figure_piece_poly.pdf")
plot(x,y ,col="gray")
lines(x, f_x, col=2, lwd=2)
abline(v=knots, lty=2)
for (i in 1:length(levels(x_cut)) ){
  subset = (x_cut==levels(x_cut)[i])
  lines(x[subset], fitted(lm(y[subset] ~ poly(x[subset],M))), lwd=2 )
}
#dev.off()

K = length(knots)
n = length(x)
B = matrix(NA,ncol=M+K+1, nrow=n)
B[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
B[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))


#pdf("Figure_smoothingmatrix.pdf")
S = B %*% solve(crossprod(B)) %*% t(B)
filled.contour(apply(S, 2, rev), asp=1, xaxt="n", yaxt="n", xlab="columns", ylab="rows", levels=c(min(S),-0.001,0.001,max(S)), col=c("red", "white", "blue"))
#dev.off()


#pdf("Figure_tpowerbasis.pdf")
plot(x,y,col="gray")
abline(v=knots, lty=2)
for (i in 1:ncol(B)) lines(x,B[,i], col=i, lwd=2)
#dev.off()

fit = lm(y ~ 0 + B)
beta_hat = coef(fit)
beta_hat
B_scaled = sapply(1:length(beta_hat),function(j) B[,j]*beta_hat[j])
#pdf("Figure_scaled_tpowerbasis.pdf")
plot(x,y ,col="gray")
abline(v=knots, lty=2)
for (i in 1:ncol(B_scaled)) lines(x,B_scaled[,i], col=i, lwd=2)
#dev.off()

#pdf("Figure_reg_spline.pdf")
plot(x,y ,col="gray")
lines(x, f_x, col=2, lwd=2)
abline(v=knots, lty=2)
y_hat = apply(B_scaled, 1, sum)
lines(x,y_hat, lwd=2)
#dev.off()



nat_spline_x <- function(x, knots){
  
  K <- length(knots)
  n <- length(x)
  xi <- knots
  
  d <- function(z, j)
  {
    out <- (x - xi[j])^3 * as.numeric(x > xi[j])
    out <- out - (x - xi[K])^3 * as.numeric(x > xi[K])
    out <- out / (xi[K] - xi[j])
    out
  }
  
  B <- matrix(0, ncol=K, nrow=n)
  B[, 1L] <- 1
  B[, 2L] <- x
  for (j in seq(1L, (K-2L)))
  {
    B[, j + 2L] <- d(x, j) - d(x, K - 1L)
  }
  B
  
}

#pdf("Figure_nat_spline.pdf")
plot(x,y ,col="gray")
lines(x, f_x, col=2, lwd=2)
abline(v=knots, lty=2)
B <- nat_spline_x(x, knots)
y_hat <- B %*% solve(crossprod(B)) %*% crossprod(B, y)
lines(x,y_hat, lwd=2)
#dev.off()


fit_ns <- lm(y ~ 0+B)
y_hat_ns <- predict(fit_ns, se = TRUE)

X = matrix(NA,ncol=M+K+1, nrow=n)
X[,1:(M+1)] = cbind(1,poly(x,M, raw=TRUE))
X[,(M+2):(M+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^M, 0))
fit_cs <- lm(y ~ 0+X)
y_hat_cs <- predict(fit_cs, se = TRUE)

#pdf("Figure_standard_error.pdf")
plot(x, y_hat_cs$se.fit, type="l", 
     ylim=c(0.02,max(y_hat_cs$se.fit)),
     ylab="Standard error", lwd=2)
lines(x, y_hat_ns$se.fit, lwd=2, col=4)
legend("top", c("cubic spline", "natural cubic spline"), lty=1, col=c(1,4), lwd=2)
#dev.off()


#pdf("Figure_smooth_spline.pdf")
overfit <- smooth.spline(x, y, all.knots=T, spar = 0)
fit_smooth <- smooth.spline(x, y, all.knots=T, cv=TRUE)
fit_smooth
plot(x,y,col="gray")
lines(x, overfit$y, col="gray")
lines(x, f_x, col=2, lwd=2)
lines(x, fit_smooth$y, col=4, lwd=2)
#dev.off()

max(svd(B)$d) / min(svd(B)$d)


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
# MCYCLE 
#---------------------------------------

rm(list=ls())

data(mcycle)
x = mcycle$times
y = mcycle$accel
