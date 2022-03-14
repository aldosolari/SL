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

D = 3
knots = seq(0.1, 0.9, by=.1)
x_cut = cut(x, c(-Inf,knots,Inf) )
#pdf("Figure_piece_poly.pdf")
plot(x,y ,col="gray")
lines(x, f_x, col=2, lwd=2)
abline(v=knots, lty=2)
for (i in 1:length(levels(x_cut)) ){
  subset = (x_cut==levels(x_cut)[i])
  lines(x[subset], fitted(lm(y[subset] ~ poly(x[subset],D))), lwd=2 )
}
#dev.off()

K = length(knots)
n = length(x)
B = matrix(NA,ncol=D+K+1, nrow=n)
B[,1:(D+1)] = cbind(1,poly(x,D, raw=TRUE))
B[,(D+2):(D+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^D, 0))
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

max(svd(B)$d) / min(svd(B)$d)

bspline_x <- function(x, knots, order=1L){
    K <- length(knots)
    k <- knots
    k <- c(rep(min(x), order + 1L), k,
           rep(max(x), order + 1L))
    b <- function(x, j, d){
      if (d == 0L) {
        return(as.numeric(k[j] <= x & x <= k[j + 1L]))
      } else {
        r <- (x - k[j]) / (k[j + d] - k[j]) * b(x, j, d - 1)
        if (is.nan(r)) r <- 0
        r <- r + (k[j + d + 1] - x) / (k[j + d + 1] - k[j + 1]) * b(x, j + 1, d - 1)
        if (is.nan(r)) r <- 0
        r
      }
    }
    B <- matrix(0, ncol=(K + order), nrow=length(x))
    for (j in seq(1L, (K + order)))
    {
      for (i in seq(1L, length(x)))
      {
        B[i, j] <- b(x[i], j, order)
      }
    }
    B
  }

B <- cbind(1,bspline_x(x, knots, order = D))
max(svd(B)$d) / min(svd(B)$d)

y_hat <- B %*% solve(crossprod(B)) %*% crossprod(B, y)

#pdf("Figure_bbasis.pdf")
plot(x,y,col="gray")
abline(v=knots, lty=2)
for (i in 1:ncol(B)) lines(x,B[,i], col=i, lwd=2)
#dev.off()

X = matrix(NA,ncol=D+K+1, nrow=n)
X[,1:(D+1)] = cbind(1,poly(x,D, raw=TRUE))
X[,(D+2):(D+K+1)] = sapply(1:K, function(k) ifelse(x >= knots[k], (x-knots[k])^D, 0))
y_t <- X %*% solve(crossprod(X)) %*% crossprod(X, y)
y_b <- B %*% solve(crossprod(B)) %*% crossprod(B, y)
max(abs(y_t - y_b))

#pdf("Figure_smooth_spline.pdf")
overfit <- smooth.spline(x, y, all.knots=T, spar = 0)
fit_smooth <- smooth.spline(x, y, all.knots=T, cv=TRUE)
plot(x,y,col="gray")
lines(x, overfit$y, col="gray")
lines(x, f_x, col=2, lwd=2)
lines(x, fit_smooth$y, col=4, lwd=2)
#dev.off()

#---------------------------------------
# MCYCLE 
#---------------------------------------

rm(list=ls())

library(MASS)
library(splines)

dataset <- mcycle %>% 
  rename(x = times, y = accel) %>%
  arrange(x)

x_grid= seq(min(dataset$x), max(dataset$x), length.out=1000)
plot(y~x, dataset, col="gray") 
fit_poly <- lm(y ~ poly(x,15), dataset)
lines(x_grid, predict(fit_poly, newdata=data.frame(x=x_grid)), lwd=2)

fit_smooth <- smooth.spline(dataset$x, dataset$y, all.knots=T)

