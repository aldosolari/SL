#=======================================
# PREDICTION, ESTIMATION AND ATTRIBUTION
#=======================================

#---------------------------------------
# Galileo data
#---------------------------------------

rm(list=ls())

time = 1:8
distance = c(33,130,298,526,824,1192,1620,2104)
D = data.frame(time, distance)

fit_Aristotle <- lm(distance ~ 0 + time, D)
fit_Galileo <- lm(distance ~ 0 + time + I(time^2), D)

#pdf("Figure_Galileo.pdf")
plot(D)
lines(fitted(fit_Aristotle))
lines(fitted(fit_Galileo), col=2)
legend("topleft", col=1:2, c("Aristotle", "Galileo"), lty=1)
#dev.off()


#---------------------------------------
# cholesterol data
#---------------------------------------

rm(list=ls())

library(readr)
library(tidyverse)

dataset <- read_table2("https://hastie.su.domains/CASI_files/DATA/cholesterol.txt") %>% 
  rename(x = compliance, y = cholesterol.decrease)
fit <- lm(y ~ poly(x,3), dataset)
x_grid <- data.frame(x=seq(min(dataset$x),max(dataset$x), length=100))
S_hat <- predict(fit, newdata = x_grid, se = TRUE)

#pdf("Figure_1.pdf")
plot(y~x, dataset, 
     pch = ".",
     xlab="normalized compliance", 
     ylab="cholesterol decrease")
lines(x_grid$x,S_hat$fit, lwd=2)
legend("topleft",
       legend=paste("Adj Rsquared = ",round(summary(fit)$adj.r.squared,3)))
ix = seq(10,90,length=11)
for (i in 1:11){
segments(x0 = x_grid$x[ix[i]], x1 = x_grid$x[ix[i]],
         y0 = S_hat$fit[ix[i]] - S_hat$se.fit[ix[i]], y1 = S_hat$fit[ix[i]] + S_hat$se.fit[ix[i]],
         col=2, lwd=2)
}
#dev.off()

summary(fit)

#---------------------------------------
# Prostate data
#---------------------------------------

rm(list=ls())

dataset <- t(read.csv("http://hastie.su.domains/CASI_files/DATA/prostmat.csv"))
dataset <- as_tibble(dataset) %>% 
  mutate(y = as.factor(c(rep("control",50),rep("cancer",52))) )

set.seed(123)
dataset_split <- rsample::initial_split(dataset, prop = 0.50, strata = y)
train <- rsample::training(dataset_split)
test  <- rsample::testing(dataset_split)

#--- random forest -------------

library(randomForest)

n_tree <- 100

set.seed(123)
rf <- randomForest(y ~ ., data = train, ntree=n_tree, importance=TRUE)

yhat <- predict(rf, test, type="vote", predict.all=TRUE)

err_by_tree = sapply(1:ncol(yhat$individual),function(i){
  apply(yhat$individual[,1:i,drop=FALSE],1,
        function(row) ifelse(mean(row=="cancer")>0.5, "cancer","control"))
})

test_errs = colMeans(err_by_tree!=test$y)
test_err = mean(ifelse(yhat$aggregate[,"cancer"] > 0.5,"cancer","control")!=test$y)

#pdf("Figure_5.pdf")
plot(1:n_tree, test_errs, type="l",
     xlab = "# of trees",
     ylab = "error rate")
legend(x=n_tree * .9, y= test_err * 2, legend = "2%", bty="n")
#dev.off()

importance_genes <- importance(rf)[,2]
important_genes <- which(importance(rf)[,2]>0)

#pdf("rf_vip_prostate.pdf")
plot(sort(importance_genes[important_genes], decreasing = T), type="h", ylab="Importance")
#dev.off()

set.seed(123)
n_rm_imp <- 100
rf_rm_imp <- randomForest(y ~ ., data = train[,-important_genes[1:n_rm_imp]], ntree=n_tree)
yhat_rm_imp <- predict(rf_rm_imp, test[,-important_genes[1:n_rm_imp]], type="vote")
mean(ifelse(yhat_rm_imp[,"cancer"] > 0.5,"cancer","control")!=test$y)


#--- gradient boosting -------------

library(xgboost)

dtrain <- xgb.DMatrix(data=as.matrix(train[,-6034]), label=(train[,6034]=="cancer")*1)
dtest <- xgb.DMatrix(data=as.matrix(test[,-6034]), label=(test[,6034]=="cancer")*1)

watchlist <- list(train=dtrain, test=dtest)

set.seed(123)
bst <- xgb.train(data=dtrain, 
                 max_depth=1, 
                 eta=0.025, 
                 nthread = 4, 
                 nrounds=n_tree, 
                 watchlist=watchlist, 
                 params = list(eval_metric  = "error"), 
                 verbose = F)

bst_test_errs <- bst$evaluation_log$test_error

#pdf("Figure_6.pdf")
plot(1:n_tree, bst_test_errs, type="l",
     xlab = "# of trees",
     ylab = "error rate")
legend(x=n_tree * .8, y= bst_test_errs[n_tree] * 1.1 , legend = "6%", bty="n")
#dev.off()

#---------------------------------------
# Prediction is Easier than Estimation
#---------------------------------------

#---------------------------------------
# The Training/Test Set Paradigm
#---------------------------------------

rm(list=ls())

dataset <- t(read.csv("http://hastie.su.domains/CASI_files/DATA/prostmat.csv"))
dataset <- as_tibble(dataset) %>% 
  mutate(y = as.factor(c(rep("control",50),rep("cancer",52))) )

train_id <- dataset[c(1:25,52:77),]
test_id <- dataset[-c(1:25,52:77),]

n_tree <- 100

set.seed(123)
rf_id <- randomForest(y ~ ., data = train_id, ntree=n_tree)

yhat_id <- predict(rf_id, test_id, type="vote", predict.all=TRUE)

err_by_tree_id = sapply(1:ncol(yhat_id$individual),function(i){
  apply(yhat_id$individual[,1:i,drop=FALSE],1,
        function(row) ifelse(mean(row=="cancer")>0.5, "cancer","control"))
})

test_errs_id = colMeans(err_by_tree_id!=test_id$y)
test_err_id = mean(ifelse(yhat_id$aggregate[,"cancer"] > 0.5,"cancer","control")!=test_id$y)

#pdf("Figure_8.pdf")
plot(1:n_tree, test_errs_id, type="l",
     xlab = "# of trees",
     ylab = "error rate")
legend(x=n_tree * .9, y= test_err_id , legend = "25%", bty="n")
#dev.off()

mean(ifelse(yhat_id$aggregate[,"cancer"] > 0.5,"cancer","control")!=test_id$y)


#---------------------------------------
# Smoothness
#---------------------------------------

rm(list=ls())

dataset <- read_table2("https://hastie.su.domains/CASI_files/DATA/cholesterol.txt") %>% 
  rename(x = compliance, y = cholesterol.decrease) %>% arrange(x)

fit_8 <- lm(y ~ poly(x,8), dataset)
train <- data.frame(y=dataset$y, x=model.matrix(fit_8)[,-1])

d <- which.min(sapply(1:8, function(d) AIC(lm(y ~ poly(x,d), dataset))))
fit <- lm(y ~ poly(x,d), dataset)

set.seed(123)
n_tree <- 200
rf <-  randomForest(y ~ ., data = train, ntree=n_tree)

#pdf("Figure_12.pdf")
plot(y~x, dataset, 
     pch=".",
     xlab = "normalized compliance",
     ylab = "cholesterol decrease",
     main = "Random Forest")
lines(dataset$x,fitted(fit), col=2, lwd=2)
lines(dataset$x,predict(rf))
#dev.off()

library(gbm)
bst <- gbm(y ~ ., data = train, distribution = "gaussian")

#pdf("Figure_12bis.pdf")
plot(y~x, dataset,
     pch=".",
     xlab = "normalized compliance",
     ylab = "cholesterol decrease", 
     main = "GBM")
lines(dataset$x,fitted(fit), col=2, lwd=2)
lines(dataset$x,fitted(fit_8), col=3)
lines(dataset$x,predict(bst))
#dev.off()

#---------------------------------------
# Traditional Methods in the Wide Data Era
#---------------------------------------

rm(list=ls())

X <- t(read.csv("http://hastie.su.domains/CASI_files/DATA/prostmat.csv"))
p = ncol(X)
pvalue <- sapply(1:ncol(X), function(i)
  t.test(x=X[1:50,i], y=X[-c(1:50),i], var.equal=T)$p.value
  )

#pdf("Figure_Bonferroni.pdf")
plot(-log10(pvalue), xlab="gene")
abline(h=-log10(0.05), lwd=2, col=3)
abline(h=-log10(0.05/p), col=2, lwd=2)
#dev.off()

library(glmnet)
y = as.factor(c(rep("control",50),rep("cancer",52)))
l <- cv.glmnet(X, y, family ="binomial")$lambda.1se
lasso <- glmnet(X, y, family ="binomial", lambda=l)

#pdf("Figure_lasso.pdf")
plot(abs(lasso$beta), xlab="gene", 
     ylab=expression(group("|",hat(beta),"|"))
     )
#dev.off()
     