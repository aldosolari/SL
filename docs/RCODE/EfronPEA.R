#---------------------------------------
# Prediction, Estimation and Attribution
#---------------------------------------

#--- Cholesterol data ------------------

rm(list=ls())

library(readr)
library(tidyverse)

dataset <- read_table2("https://hastie.su.domains/CASI_files/DATA/cholesterol.txt") %>% 
  rename(x = compliance, y = cholesterol.decrease)
fit <- lm(y ~ poly(x,3), dataset)
x_grid <- data.frame(x=seq(min(dataset$x),max(dataset$x), length=100))
Shat <- predict(fit, newdata = x_grid, se = TRUE)
plot(y~x, dataset, 
     pch = ".",
     xlab="normalized compliance", 
     ylab="cholesterol decrease")
lines(x_grid$x,Shat$fit, lwd=2)
legend("topleft",
       legend=paste("Adj Rsquared = ",round(summary(fit)$adj.r.squared,3)))
ix = seq(10,90,length=11)
for (i in 1:11){
segments(x0 = x_grid$x[ix[i]], x1 = x_grid$x[ix[i]],
         y0 = Shat$fit[ix[i]] - Shat$se.fit[ix[i]], y1 = Shat$fit[ix[i]] + Shat$se.fit[ix[i]],
         col=2, lwd=2)
}
