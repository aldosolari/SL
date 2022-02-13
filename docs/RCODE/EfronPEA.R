#=======================================
# PREDICTION, ESTIMATION AND ATTRIBUTION
#=======================================

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
library(vip)

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

#pdf("rf_vip.pdf")
rf_fit_best %>% 
  extract_fit_parsnip() %>% 
  vip() 
#dev.off()

rf_res_best <- rf_wflow_best %>% 
  fit_resamples(resamples = dataset_folds, control = keep_pred) 

compare_models <- as_workflow_set(rf = rf_res_best) %>% 
  bind_rows(compare_models)

collect_metrics(compare_models) %>% 
  mutate(wflow_id = gsub("(workflow_variables_)", "", wflow_id)) %>%
  filter(.metric == "accuracy")

rf_accuracy <- collect_metrics(rf_res_best ) %>% filter(.metric == "accuracy") %>% select(mean)
  
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

#pdf("rf_vip_prostate.pdf")
varImpPlot(rf, type=1)
#dev.off()

#--- gradient boosting -------------

library(xgboost)

dtrain <- xgb.DMatrix(data=as.matrix(train[,-6034]), label=(train[,6034]=="cancer")*1)
dtest <- xgb.DMatrix(data=as.matrix(test[,-6034]), label=(test[,6034]=="cancer")*1)

watchlist <- list(train=dtrain, test=dtest)

set.seed(123)
bst <- xgb.train(data=dtrain, 
                 max_depth=4, 
                 eta=0.01, 
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
legend(x=n_tree * .8, y= bst_test_errs[n_tree] * 1.1 , legend = "4%", bty="n")
#dev.off()
