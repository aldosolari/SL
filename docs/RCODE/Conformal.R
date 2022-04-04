#=======================================
# CONFORMAL PREDICTION
#=======================================

rm(list=ls())

alpha = 0.1
n = 100

set.seed(123)
x = sort(runif(n,-5,5))
y = 1/4 * (x+4) * (x+1) * (x-2) + rnorm(n, mean = 1, sd = 2)
train <- data.frame(x,y)

x_new = runif(1,-5,5)
C = predict(lm(y ~ poly(x,degree=3), train), 
            newdata=data.frame(x=x_new), 
            interval = "prediction",
            level = 1-alpha)
y_new = 1/4 * (x_new+4) * (x_new+1) * (x_new-2) + rnorm(1, mean = 1, sd = 2)


#pdf("Figure_prediction_interval.pdf")
plot(y~x,train)
lines(x, 1/4 * (x+4) * (x+1) * (x-2), lwd=2)
rug(x_new)
segments(x0=x_new, x1=x_new, y0=C[,2], y1=C[,3], col=2, lwd=2)
points(x_new,y_new,pch=19)
#dev.off()

B <- 1000
coverage <- vector()
for (i in 1:B){
  x_train = runif(n,-5,5)
  y_train = 1/4 * (x_train+4) * (x_train+1) * (x_train-2) + rnorm(n, mean = 1, sd = 2)
  C = predict(lm(y ~ poly(x,degree=3), data.frame(x_train,y_train)), 
              newdata=data.frame(x=x_new), 
              interval = "prediction",
              level = 1-alpha)
  y_new = 1/4 * (x_new+4) * (x_new+1) * (x_new-2) + rnorm(1, mean = 1, sd = 2)
coverage[i] = C[,2] <= y_new & y_new <= C[,3]
}
mean(coverage)

#---------------------------------------
# MODEL MISS-SPECIFICATION
#---------------------------------------

C = predict(lm(y ~ x, train), 
            newdata=data.frame(x=train$x), 
            interval = "prediction",
            level = 1-alpha)

#pdf("Figure_wrong_specification.pdf")
plot(y~x,train)
lines(x, 1/4 * (x+4) * (x+1) * (x-2), lwd=2)
polygon(c(x,rev(x)), 
        c(C[,2],rev(C[,3])),
        col=rgb(1, 0, 0,0.5), border=NA)
#dev.off()

B <- 1000
coverage_mat <- matrix(NA,nrow=B,ncol=n)
x_mat <- coverage_mat
for (i in 1:B){
  x_train = runif(n,-5,5)
  y_train = 1/4 * (x_train+4) * (x_train+1) * (x_train-2) + rnorm(n, mean = 1, sd = 2)
  x_new = runif(n,-5,5)
  x_mat[i,] <- x_new
  C = predict(lm(y ~ x, data.frame(x=x_train,y=y_train)), 
              newdata=data.frame(x=x_new), 
              interval = "prediction",
              level = 1-alpha)
  y_new = 1/4 * (x_new+4) * (x_new+1) * (x_new-2) + rnorm(1, mean = 1, sd = 2)
  coverage_mat[i,] = C[,2] <= y_new & y_new <= C[,3]
}

mean(coverage_mat)
coverage_tab <- aggregate(c(coverage_mat), by=list(cut(x_mat,50)), mean)

#pdf("Figure_coverage_wrong_specification.pdf")
barplot(coverage_tab[,2], names.arg=coverage_tab[,1], ylim=c(0,1), ylab="Coverage", xlab="x", xaxt="n")
abline(h=1-alpha, lwd=2,col=2)
#dev.off()


#---------------------------------------
# SPLIT CONFORMAL
#---------------------------------------

split_conformal = function(x, y, x_new, m, alpha=0.1,
                           split=NULL, seed=NULL){
  
  require(randomForest)
  
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x_new = matrix(x_new,ncol=p)
  n_new = nrow(x_new)
  
  if (!is.null(split)) I = split
  else {
    if (!is.null(seed)) set.seed(seed)
    I = sample(1:n,m)
  }
  L = (1:n)[-I]

  fit = randomForest(x=x[L,,drop=F],y=y[L])
  y_new = matrix(predict(fit,x_new),nrow=n_new)
  
  res =  abs(y[I] - predict(fit,x[I,,drop=F]))
  o = order(res)
  c = ceiling((1-alpha)*(m+1))
  r = res[o][c]
  
  lo = up = vector()
  for (i in 1:n_new) {
    lo[i] = y_new[i] - r 
    up[i] = y_new[i] + r 
  }
  
  return(list(lo=lo,up=up))
}

x_new = seq(-5,5, length.out = 1000)
C = split_conformal(x, y, x_new,
                    alpha = 0.1,
                    m = 49)

#pdf("Figure_random_forest.pdf")
plot(y~x,train)
lines(train$x, 1/4 * (train$x+4) * (train$x+1) * (train$x-2), lwd=2)
polygon(c(x_new,rev(x_new)), 
        c(C$lo,rev(C$up)),
        col=rgb(1, 0, 0,0.5), border=NA)
#dev.off()

B <- 1000
coverage_mat <- matrix(NA,nrow=B,ncol=n)
x_mat <- coverage_mat

for (i in 1:B){
  x = runif(n,-5,5)
  y = 1/4 * (x+4) * (x+1) * (x-2) + rnorm(n, mean = 1, sd = 2)
  x_new = runif(n,-5,5)
  x_mat[i,] <- x_new
  C = split_conformal(x, y, x_new,
                      alpha = 0.1,
                      m = 49)
  y_new = 1/4 * (x_new+4) * (x_new+1) * (x_new-2) + rnorm(1, mean = 1, sd = 2)
  coverage_mat[i,] = C$lo <= y_new & y_new <= C$up
}

mean(coverage_mat)
coverage_tab <- aggregate(c(coverage_mat), by=list(cut(x_mat,50)), mean)

#pdf("Figure_coverage_random_forest.pdf")
barplot(coverage_tab[,2], names.arg=coverage_tab[,1], ylim=c(0,1), ylab="Coverage", xlab="x", xaxt="n")
abline(h=1-alpha, lwd=2,col=2)
#dev.off()

#---------------------------------------
# ORACLE
#---------------------------------------

mu_x = 1
sigma_x = 1
mu_y = 2
sigma_y = 1
rho = 0.8

set.seed(123)
n = 10^3
x_i = sort( rnorm(n, mean=mu_x, sd=sigma_x) )
mu_yIx = mu_y + rho * (sigma_x / sigma_y) * (x_i - mu_x)
sigma_yIx = sqrt( (sigma_y)^2 * (1-rho^2) )
y_i = rnorm(n, mu_yIx, sd = sigma_yIx)

q1_x_i = qnorm(alpha/2, mean=mu_yIx, sd = sigma_yIx )
q2_x_i = qnorm(1-alpha/2, mean=mu_yIx, sd = sigma_yIx )

#pdf("Figure_oracle.pdf")
plot(x_i,y_i, xlab="x", ylab="y")
polygon(c(x_i,rev(x_i)), 
        c(q1_x_i,rev(q2_x_i)),
        col=rgb(1, 0, 0,0.5), border=NA)
#dev.off()

coverage = y_i >= q1_x_i & y_i <= q2_x_i
coverage_tab <- aggregate(coverage, by=list(cut(x_i,quantile(x_i, probs=seq(0,1,0.1)))), mean)
barplot(coverage_tab[,2], names.arg=coverage_tab[,1], ylim=c(0,1), ylab="Coverage", xlab="x", xaxt="n")


#---------------------------------------
# QUANTILE CONFORMAL
#---------------------------------------

split_conformal_quantile = function(x, y, x_new, m,
                    alpha=0.1, 
                    gamma = alpha/2, 
                    split=NULL, seed=NULL) {
  
  require(quantregForest)
  
  x = as.matrix(x)
  y = as.numeric(y)
  n = nrow(x)
  p = ncol(x)
  x_new = matrix(x_new,ncol=p)
  n_new = nrow(x_new)
  
  if (!is.null(split)) I = split
  else {
    if (!is.null(seed)) set.seed(seed)
    I = sample(1:n,m)
  }
  L = (1:n)[-I]
  n_L = length(L)
  n_I = length(I)
  
  fit = quantregForest(x=x[L,,drop=F],y=y[L], nthreads=16)
  y_new = matrix(predict(fit,x_new, what=c(gamma,1-gamma)),nrow=n_new)
  
  res = apply( cbind(y[I],-y[I]) + matrix(predict(fit,x[I,,drop=F], what=c(1-gamma,gamma)),nrow=n_I) %*% diag(c(-1,1)),1,max )
  
  o = order(res)
  c = ceiling((1-alpha)*(m+1))
  r = res[o][c]
  
  lo = up = vector()
  
  for (i in 1:n_new) {
    lo[i] = y_new[i,1] - r 
    up[i] = y_new[i,2] + r 
  }
  
  return(list(lo=lo,up=up))
}


set.seed(123)

n = 100
x = sort(runif(n,0,2*pi))
y = sin(x) + x*pi/30*rnorm(n)

x_new = seq(0,2*pi,length=1000)

C = split_conformal_quantile(x, y, x_new,
                             alpha = 0.1,
                             m = 49)

#pdf("Figure_conformal_quantile.pdf")
plot(y~x)
lines(x, sin(x), lwd=2)
polygon(c(x_new,rev(x_new)), 
        c(C$lo,rev(C$up)),
        col=rgb(1, 0, 0,0.5), border=NA)
#dev.off()

#---------------------------------------
# MULTI SPLIT
#---------------------------------------

multi_split <- function(C_mat, tau=0.5){

 n0 = nrow(C_mat)/2
 B = ncol(C_mat)
 Y = cbind(C_mat[1:n0,,drop=F],C_mat[(n0+1):(2*n0),,drop=F])
 H = matrix(rep(rep(1:0,each=B),n0), byrow=T, ncol=2*B)
 tr <- tau*B + .001
 lo <- rep(NA,n0)
 up <- rep(NA,n0)  
  
for (i in 1:n0){
  
  y = Y[i,,drop=FALSE]
  h = H[i,,drop=FALSE]
  o = order(y,2-h)
  ys <- y[o]
  hs <- h[o]
  
  count <- 0
  leftend <- 0
  
  for (j in 1:(2*B) ){
    if ( hs[j]==1 ) {
      count <- count + 1
      if ( count > tr && (count - 1) <= tr) { 
        leftend <- ys[j] 
      }
    } else {
      if ( count > tr && (count - 1) <= tr) {
        rightend <- ys[j]
        lo[i] <- leftend
        up[i] <- rightend
      }
      count <- count - 1
    }
  }
}  

 return(list(lo=lo,up=up))
}


set.seed(123)

n = 100
x = sort(runif(n,0,2*pi))
y = sin(x) + x*pi/30*rnorm(n)

n_new = 1
x_new = runif(n_new,0,2*pi)

B = 10
C_mat = replicate(B,unlist(split_conformal_quantile(x, y, x_new,
                                                    alpha = 0.1,
                                                    m = 49)))

C_multi = multi_split(C_mat,tau=0.5)

plot(NULL, xlim=c(1,B), ylim=c(min(C_mat[1,]), max(C_mat[2,])), ylab="C", xlab="split")
for (b in 1:B){
  segments(x0=b,x1=b,y0=C_mat[1,b],y1=C_mat[2,b])
}
abline(h=C_multi$lo)
abline(h=C_multi$up)

