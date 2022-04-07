#=======================================
# CLASSICAL VS HIGHDIMENSIONAL THEORY
#=======================================


#---------------------------------------
# Logistic regression
#---------------------------------------

rm(list=ls())

n <- 4000
p <- 800

beta = c(rep(10,100),rep(-10,100),rep(0,600))
set.seed(123)
X = matrix(rnorm(n*p, mean=0, sd = 1/sqrt(n)), ncol=p)
lp = X%*%beta
pr_1 <- exp(lp)/(1+exp(lp))
y <- rbinom(n,size=1,prob=pr_1)

fit <- glm(y~X,family="binomial")

#pdf("Figure_logistic_coef.pdf")
plot(1:p,beta, type="p", col=2, xlab="index", ylab="coefficients (true and fitted)", ylim=c(-30,30))
points(1:p,fit$coefficients[-1])
#dev.off()

set.seed(123)
beta = rnorm(p,3,sqrt(16))

lp = X%*%beta
pr_1 <- exp(lp)/(1+exp(lp))
y <- rbinom(n,size=1,prob=pr_1)

fit <- glm(y~X,family="binomial")

#pdf("Figure_logistic_slope.pdf")
plot(beta,fit$coefficients[-1],  xlab="True coefficient", ylab="Estimated coefficient", asp=1)
abline(a=0,b=1,lwd=2,col=2)
#dev.off()

set.seed(123)
X_new <- matrix(rnorm(n*p, mean=0, sd = 1/sqrt(n)), ncol=p)

lp_true <- X_new %*% beta
p_1_true <- exp(lp_true)/(1+exp(lp_true))
p_1_hat <- predict(fit, newdata=data.frame(X), "response")
ix <- sort(p_1_true, index.return=TRUE)$ix

#pdf("Figure_logistic_prob.pdf")
plot(p_1_true[ix], type="l", col=2, lwd=2, ylab="Probability (true and predicted)", xaxt="n", xlab="")
points(p_1_hat[ix], pch=".")
#dev.off()


set.seed(123)
beta = c(rep(0,p/2), rnorm(p/2,7,1))

lp = X%*%beta
pr_1 <- exp(lp)/(1+exp(lp))
y <- rbinom(n,size=1,prob=pr_1)

fit <- glm(y~X,family="binomial")

#pdf("Figure_logistic_pvalues.pdf")
hist(summary(fit)$coeff[2:((p/2)+1),4], freq=F, 20, xlab="P-values", main="")
#dev.off()


#---------------------------------------
# Linear discriminant analysis
#---------------------------------------

rm(list=ls())

library(MASS)

n <- 800
p <- 400

sim_err_hat = function(gamma, p, delta){
  mu_A = rep(0,p)
  mu_B = rep(gamma/sqrt(p),p)
  n = p / delta
  x_A = mvrnorm(n,mu_A, Sigma = diag(rep(1,p)))
  x_B = mvrnorm(n,mu_B, Sigma = diag(rep(1,p)))
  mu_A_hat = colMeans(x_A)
  mu_B_hat = colMeans(x_B)
  SD = sqrt( crossprod( mu_A_hat - mu_B_hat ) )
  pr_A = pnorm(0, crossprod(mu_A_hat - mu_B_hat, mu_A - (mu_A_hat + mu_B_hat)/2), sd = SD )
  pr_B = 1-pnorm(0, crossprod(mu_A_hat - mu_B_hat, mu_B - (mu_A_hat + mu_B_hat)/2), sd = SD )
  Err_hat = pr_A/2 + pr_B/2
  return(Err_hat)
}


gammas = seq(1,2,length.out = 6)
delta = p/n

Err = pnorm(-gammas/2)
Err_hd = pnorm(-gammas^2/(2*sqrt(gammas^2 + 2*delta)))
plot(gammas, Err_hd, type="l", ylab="Error probability",
     ylim=c(0.1,0.4), xlab=expression(gamma))
lines(gammas, Err, lty=2)
legend(legend=c("hd","classical"), "topright", lty=1:2, bty="n")
Err_hat = vector()

set.seed(123)
B = 1
for (i in 1:length(gammas)){
  Err_hat[i] = mean(replicate(B, sim_err_hat(gamma=gammas[i],p=p, delta=p/n)))
}
points(gammas, Err_hat)


deltas = c(0.05,.2,.5,.8,1)

Err = pnorm(-1)
Err_hd = pnorm(-2^2/(2*sqrt(2^2 + 2*deltas)))
plot(deltas, Err_hd, type="l", ylab="Error probability",
     ylim=c(0.154,0.22),xlab=expression(delta))
abline(h=Err, lty=2)
Err_hat = vector()

set.seed(123)
B = 1
for (i in 1:length(deltas)){
  Err_hat[i] = mean(replicate(B, sim_err_hat(gamma=2,delta=deltas[i],p=p)))
}
points(deltas, Err_hat)



sim_err_hat_s = function(gamma, p, s, delta){
  mu = rep(0,p)
  mu[1:s]<- gamma/sqrt(s)
  mu_A = mu/2
  mu_B = -mu/2
  n = p / delta
  x_A = mvrnorm(n,mu_A, Sigma = diag(rep(1,p)))
  x_B = mvrnorm(n,mu_B, Sigma = diag(rep(1,p)))
  lambda <- sqrt(2*log(p)/n)
  mu_A_hat = colMeans(x_A)
  mu_B_hat = colMeans(x_B)
  mu_A_hat[abs(mu_A_hat) < lambda]<-0
  mu_B_hat[abs(mu_B_hat) < lambda]<-0
  SD = sqrt( t(mu_A_hat - mu_B_hat) %*% (mu_A_hat - mu_B_hat) )
  pr_A = pnorm(0, crossprod(mu_A_hat - mu_B_hat, mu_A - (mu_A_hat + mu_B_hat)/2), sd = SD )
  pr_B = 1-pnorm(0, crossprod(mu_A_hat - mu_B_hat, mu_B - (mu_A_hat + mu_B_hat)/2), sd = SD )
  Err_hat = pr_A/2 + pr_B/2
  return(Err_hat)
}


Err = pnorm(-gammas/2)
Err_hd = pnorm(-gammas^2/(2*sqrt(gammas^2 + 2*delta)))
plot(gammas, Err_hd, type="l", ylab="Error probability",
     ylim=c(0.1,0.4), xlab=expression(gamma))
lines(gammas, Err, lty=2)
legend(legend=c("hd","classical"), "topright", lty=1:2, bty="n")
Err_hat = vector()

set.seed(123)
B = 1
for (i in 1:length(gammas)){
  Err_hat[i] = mean(replicate(B, sim_err_hat_s(gamma=gammas[i],p=p, s=5, delta=p/n)))
}
points(gammas, Err_hat)




