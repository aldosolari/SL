#=======================================
# JAMES-STEIN ESTIMATION
#=======================================

#---------------------------------------
# Upper bound Risk JS
#---------------------------------------

rm(list=ls())

p <- 10
mu_ss = seq(0,2*p, length=100)
bound_JS <- p - (p-2) / (1 + mu_ss / (p-2) )
oracle <- p*mu_ss / (p + mu_ss)

#pdf("Figure_JSbound.pdf")
plot(mu_ss, bound_JS, type="l", 
     ylim = c(0,p),
     xlab = expression("||"*mu*"||" ^2),
     ylab = "Risk", 
     main = paste("p = ", p))
abline(h=p, lty=2)
lines(mu_ss, oracle, lty=3)
legend("bottomright", c("MLE","JS (bound)","oracle"), lty=c(2,1,3))
#dev.off()


#---------------------------------------
# Simulation: Risk
#---------------------------------------

rm(list=ls())

library(xtable)

sim_risk <- function(p){
  q = round(p/2)
  mu = c(rep(sqrt(p/q),q), rep(0,p-q))
  x <- rnorm(p, mean = mu )
  mu_hat <- x
  mu_hat_JS <- ( 1- ( (p-2)/sum(x^2) ) )*x
  loss <- (mu - mu_hat)^2 
  loss_JS <- (mu - mu_hat_JS)^2
  return(c(MLE=loss,JS=loss_JS))
}

B <- 20000
p = 5
set.seed(123)
res = replicate(B, sim_risk(p))

MLE = c(mean(colSums(res[grepl("MLE", rownames(res)),])),
       apply(res[grepl("MLE", rownames(res)), ],1,mean))
JS = c(mean(colSums(res[grepl("JS", rownames(res)),])),
        apply(res[grepl("JS", rownames(res)), ],1,mean))

table_res = rbind(MLE,JS)
colnames(table_res) <- c("Risk",paste0("Risk",1:p))
xtable(table_res)

#---------------------------------------
# Simulation: Bayes Risk
#---------------------------------------

rm(list=ls())

sim_Bayes_risk <- function(p, tau2){
  mu = rnorm(p, mean=0, sd=sqrt(tau2))
  x <- rnorm(p, mean = mu )
  mu_hat <- x
  mu_hat_B <- (tau2/(1+tau2))*x
  mu_hat_JS <- ( 1- ( (p-2)/sum(x^2) ) )*x
  loss <- (mu - mu_hat)^2 
  loss_B <- (mu - mu_hat_B)^2
  loss_JS <- (mu - mu_hat_JS)^2
  return(c(MLE=loss,BAYES= loss_B, JS=loss_JS))
}

B <- 20000
p = 5
tau2 = 2
set.seed(123)
res = replicate(B, sim_Bayes_risk(p, tau2))

MLE = c(mean(colSums(res[grepl("MLE", rownames(res)),])),
        apply(res[grepl("MLE", rownames(res)), ],1,mean))
BAYES = c(mean(colSums(res[grepl("BAYES", rownames(res)),])),
       apply(res[grepl("BAYES", rownames(res)), ],1,mean))
JS = c(mean(colSums(res[grepl("JS", rownames(res)),])),
       apply(res[grepl("JS", rownames(res)), ],1,mean))

table_res = rbind(MLE,BAYES,JS)
colnames(table_res) <- c("Bayes Risk",paste0("B.Risk",1:p))
xtable(table_res)

#---------------------------------------
# Baseball data
#---------------------------------------

rm(list=ls())

dataset <- read_table2("https://hastie.su.domains/CASI_files/DATA/baseball.txt")[,-1]
xtable(dataset)

y <- dataset$MLE
p <- length(y)
n <- 90
x <- 2*sqrt(n+0.5)*asin( sqrt((n*y + 0.375)/(n+0.75)) )

xbar <- mean(x)
S <- sum((x-xbar)^2)
mu_hat_JS <- xbar + (1-((p-3)/S))*(x - xbar)

y_JS <- (1/n) *( (n+0.75) * (sin(mu_hat_JS/(2*sqrt(n+0.5)) ) )^2 - 0.375 )
y_TRUE <- dataset$TRUTH

sum((y-y_TRUE)^2)
sum((y-y_JS)^2)

sum( (y-y_TRUE)^2 > (y-y_JS)^2 )

#pdf("Figure_baseball.pdf")
plot(c(y[1],y_JS[1],y_TRUE[1]),2:0, type="b", lty=2, ylim=c(0,2), xlim=range(y), yaxt="n", ylab="", xlab="Batting average")
for (i in 2:p) lines(c(y[i],y_JS[i],y_TRUE[i]),2:0, lty=2, type="b")
axis(2, at=0:2, labels=c("TRUE","JS","MLE"))
#dev.off()
