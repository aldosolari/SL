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
Shat <- predict(fit, newdata = x_grid, se = TRUE)

#pdf("Figure_1.pdf")
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
#dev.off()

summary(fit)

#---------------------------------------
# birthwt data
#---------------------------------------

rm(list=ls())

library(MASS)
library(tidymodels)
library(rpart.plot)

#--- pre-processing of Venables & Ripley -----

dataset <- birthwt[,-10]
dataset$race <- factor(dataset$race, labels = c("white", "black", "other"))
dataset$ptd <- (dataset$ptl > 0)*1
dataset$ftv <- factor(dataset$ftv)
levels(dataset$ftv)[-(1:2)] <- "2+"
dataset$low <- factor(dataset$low)
dataset$smoke <- (dataset$smoke > 0)*1
dataset$ht <- (dataset$ht > 0)*1
dataset$ui <- (dataset$ui > 0)*1

#--- logistic regression -------------

lr_model <- logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") 

lr_wflow <- workflow() %>%
  add_model(lr_model) %>%
  add_formula(low ~ .)

lr_fit <- fit(lr_wflow, dataset)

tidy(lr_fit)

augment(lr_fit, new_data = dataset) %>%
  conf_mat(truth = low, estimate = .pred_class)

augment(lr_fit, new_data = dataset) %>%
  accuracy(truth = low, estimate = .pred_class)

set.seed(123)
dataset_folds <- vfold_cv(dataset, v = 10)

keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)

lr_res <- lr_wflow %>% 
  fit_resamples(resamples = dataset_folds, control = keep_pred) 

lr_res %>% collect_metrics() %>%
  filter(.metric == "accuracy")

#--- classification tree -------------

tree_model <- decision_tree() %>%
  set_engine("rpart") %>%
  set_mode("classification") 

tree_wflow <- workflow() %>%
  add_model(tree_model) %>%
  add_formula(low ~ .)

tree_fit <- fit(tree_wflow, dataset)

#pdf("Figure_3.pdf")
tree_fit %>%
  extract_fit_engine() %>%
  rpart.plot(roundint=FALSE)
#dev.off()

augment(tree_fit, new_data = dataset) %>%
  accuracy(truth = low, estimate = .pred_class)

tree_res <- tree_wflow %>% 
  fit_resamples(resamples = dataset_folds, control = keep_pred) 

compare_models <- 
  as_workflow_set(lr = lr_res, tree = tree_res) 

collect_metrics(compare_models) %>% 
  mutate(wflow_id = gsub("(workflow_variables_)", "", wflow_id)) %>%
  filter(.metric == "accuracy")

#--- random forest -------------

n_trees <- 1000

rf_model <- rand_forest(
  trees = n_trees, 
  mtry=tune()
  ) %>%
  set_engine("randomForest", importance = TRUE) %>%
  set_mode("classification") 

rf_wflow <- workflow() %>%
  add_model(rf_model) %>%
  add_formula(low ~ .) 

set.seed(123)
rf_tuned <- 
  rf_wflow %>% 
  tune_grid(
    resamples = dataset_folds,
    grid = 3
  )

rf_tuned %>% 
  collect_metrics() %>% 
  filter(.metric == "accuracy")

rf_best <- rf_tuned %>%
  select_best("accuracy")

rf_wflow_best <- 
  rf_wflow %>% 
  finalize_workflow(rf_best)

rf_fit_best <- fit(rf_wflow_best, dataset)

augment(rf_fit_best, new_data = dataset) %>%
  accuracy(truth = low, estimate = .pred_class)

rf_res_best <- rf_wflow_best %>% 
  fit_resamples(resamples = dataset_folds, control = keep_pred) 

compare_models <- as_workflow_set(rf = rf_res_best) %>% 
  bind_rows(compare_models)

collect_metrics(compare_models) %>% 
  mutate(wflow_id = gsub("(workflow_variables_)", "", wflow_id)) %>%
  filter(.metric == "accuracy")

rf_accuracy <- collect_metrics(rf_res_best) %>% filter(.metric == "accuracy") %>% dplyr::select(mean)
  
#pdf("Figure_4.pdf")
OOB <- rf_fit_best[["fit"]][["fit"]][["fit"]][["err.rate"]][,"OOB"]
plot(OOB, type="l", xlab="# of trees", ylab="prediction error")
abline(h=1-rf_accuracy, lty="dashed")
legend(x = n_trees * 0.9, y=1-rf_accuracy, legend="31%", bty="n")
#dev.off()

#---------------------------------------
# Prostate data
#---------------------------------------

rm(list=ls())

dataset <- t(read.csv("http://hastie.su.domains/CASI_files/DATA/prostmat.csv"))
dataset <- as_tibble(dataset) %>% 
  mutate(y = as.factor(c(rep("control",50),rep("cancer",52))) )

set.seed(123)
dataset_split <- initial_split(dataset, prop = 0.50, strata = y)
train <- training(dataset_split)
test  <- testing(dataset_split)

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

rm(list=ls())

smallsim <- function(mu){
  x <- rnorm(25, mu)
  xbar <- mean(x)
  xtilde <- median(x)
  ee_xbar <- (xbar - mu)^2
  ee_xtilde <- (xtilde - mu)^2
  X <- rnorm(1, mu)
  pe_xbar <- (xbar - X)^2
  pe_xtilde <- (xtilde - X)^2
  return(c(ee_xbar,ee_xtilde,pe_xbar,pe_xtilde))
}

set.seed(123)
B <- 1000
res <- replicate(B,smallsim(1))
mean( res[2,] ) / mean( res[1,] )
mean( res[4,] ) / mean( res[3,] )

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
# Concept drift
#---------------------------------------

rm(list=ls())

n <- 400
p <- 200
set.seed(123)
X <- matrix(rnorm(n*p), ncol=p)
y = rep(c("treatment","control"), n/2)

Z <- matrix(0, nrow=n,ncol=p)
for (j in 1:p){
  episode <- which(rbinom(n,size=1,prob=1/n)==1)
  episode_days <- unique(unlist(sapply(episode, function(x) x:min((x+30),n)) ))
  Z[episode_days,j] <- 1
  if (length(episode_days)>0){
    segno = sample(c(-1,1), size=1)
    X[episode_days[episode_days%%2==1],j] <- X[episode_days[episode_days%%2==1],j] + segno*2
    X[episode_days[episode_days%%2==0],j] <- X[episode_days[episode_days%%2==0],j] - segno*2
  }
}

#pdf("Figure_10.pdf")
image(Z, col=c(0,1), asp=p/n, xlab="days (subjects)", ylab="genes", xaxt="n", yaxt="n")
abline(v=0.8)
#dev.off()

dataset <- data.frame(y=as.factor(y),X)
n_tree <- 200

#--- random forest random split ----------

id_random <- sample(1:n,size=0.8*n)
train_random <- dataset[id_random,]
test_random <- dataset[-id_random,]

rf_random <- randomForest(y ~ ., data = train_random, ntree=n_tree)

yhat_random <- predict(rf_random, test_random, type="vote", predict.all=TRUE)

err_by_tree_random = sapply(1:ncol(yhat_random$individual),function(i){
  apply(yhat_random$individual[,1:i,drop=FALSE],1,
        function(row) ifelse(mean(row=="treatment")>0.5, "treatment","control"))
})

test_errs_random = colMeans(err_by_tree_random!=test_random$y)
test_err_random = mean(ifelse(yhat_random$aggregate[,"treatment"] > 0.5,"treatment","control")!=test_random$y)

#pdf("Figure_11.pdf")
plot(1:n_tree, test_errs_random, type="l",
     ylim = c(0,.5),
     xlab = "# of trees",
     ylab = "error rate")
legend(x=n_tree * .8, y= test_err_random, legend = "6%", bty="n")
#dev.off()

#--- random forest ordered split ----------

id_ordered <- 1:(n*0.8)
train_ordered <- dataset[id_ordered,]
test_ordered <- dataset[-id_ordered,]

rf_ordered <- randomForest(y ~ ., data = train_ordered, ntree=n_tree)

yhat_ordered <- predict(rf_ordered, test_ordered, type="vote", predict.all=TRUE)

err_by_tree_ordered = sapply(1:ncol(yhat_ordered$individual),function(i){
  apply(yhat_ordered$individual[,1:i,drop=FALSE],1,
        function(row) ifelse(mean(row=="treatment")>0.5, "treatment","control"))
})

test_errs_ordered = colMeans(err_by_tree_ordered!=test_ordered$y)
test_err_ordered = mean(ifelse(yhat_ordered$aggregate[,"treatment"] > 0.5,"treatment","control")!=test_ordered$y)

#pdf("Figure_11bis.pdf")
plot(1:n_tree, test_errs_ordered, type="l",
     ylim = c(0,.5),
     xlab = "# of trees",
     ylab = "error rate")
legend(x=n_tree * .8, y= test_err_ordered, legend = "17%", bty="n")
#dev.off()

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
     