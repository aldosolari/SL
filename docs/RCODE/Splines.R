#=======================================
# SPLINES
#=======================================


#---------------------------------------
# DOPPLER FUNCTION 
#---------------------------------------



#---------------------------------------
# EXAMPLE
#---------------------------------------

rm(list=ls())




library(MASS)

dataset <- mcycle %>% 
  rename(x = times, y = accel) %>%
  arrange(x)

x_grid= seq(min(dataset$x), max(dataset$x), length.out=1000)
d = 15
plot(y~x, dataset, col="gray") 
fit <- lm(y ~ poly(x,d), dataset)
lines(x_grid, predict(fit, newdata=data.frame(x=x_grid)), lwd=2)


n <- 500
sigma <- 0.3
d <- 15
  
set.seed(123)
x = sort(runif(n)) 
f_x = sin(2*(4*x-2)) + 2*exp(-(16^2)*((x-.5)^2))
y = f_x + rnorm(n, mean=0, sd=sigma)
plot(x, y, col="gray") 
lines(x, f_x, col=2, lwd=2)
lines(x, fitted(lm(y ~ poly(x,d))), lwd=2)