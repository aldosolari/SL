#=======================================
# SPARSITY
#=======================================


#---------------------------------------
# TOY EXAMPLE
#---------------------------------------

rm(list=ls())

y = seq(-8,8,length.out = 401)
lambda = 1
mu_hat_0 = y*(abs(y) > sqrt(2*lambda))
mu_hat_1 = (y+sign(lambda-y)*lambda)*(abs(y) > lambda)
mu_hat_2 = (1/(1+2*lambda))*y

#pdf("Figure_toy.pdf")
plot(y, mu_hat_0, pch=19, asp=1, col=2, xlim=c(-4,4), ylim=c(-4,4), ylab=expression(hat(mu)))
abline(a=0,b=1, lty=3)
abline(h=0, lty=3)
abline(v=0, lty=3)
points(y, mu_hat_1, pch=19, col=3 )
points(y, mu_hat_2, pch=19, col=4 )
legend("topleft", c("l0","l1","l2"), col=c(2,3,4), pch=19)
#dev.off()

#---------------------------------------
# ORTHOGONAL CASE
#---------------------------------------

rm(list=ls())

n = 400
p = 4
set.seed(123)
Z <- matrix(rnorm(n*p), ncol = p)
X <- svd(Z)$u
round(crossprod(X), 10)
beta = c(2,rep(0,p-1))
y <- X %*% beta + rnorm(n)
Xty <- crossprod(X, y)
lambdas = seq(0,7, length.out = 500)
beta_hat_0 = sapply(1:length(lambdas), function(i)
             Xty * (abs(Xty) > sqrt(2*lambdas[i])) )
beta_hat_1 = sapply(1:length(lambdas), function(i)
  (Xty+sign(lambdas[i]-Xty)*lambdas[i])*(abs(Xty) > lambdas[i]) )
beta_hat_2 = sapply(1:length(lambdas), function(i)
             Xty * (1/(1+2*lambdas[i])) )

#pdf("Figure_orthogonal.pdf")
matplot(lambdas, t(beta_hat_2), type="l", lty=1, lwd=2, col=4,
        ylab=expression(hat(beta)), xlab=expression(lambda))
matlines(lambdas, t(beta_hat_1), col=3, lty=1,lwd=2)
matlines(lambdas, t(beta_hat_0), col=2, lty=1, lwd=2)
points(rep(0,p),Xty, pch=19)
legend("topright", c("BSS","Ridge","Lasso"), col=c(2,4,3), lty=1, lwd=2)
#dev.off()

#---------------------------------------
# PROSTATE
#---------------------------------------

rm(list=ls())

library(readr)
library(tidyverse)

dataset <- read_delim("https://hastie.su.domains/ElemStatLearn/datasets/prostate.data", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

train <- dataset %>%  filter(train) %>% select(-X1,-train) %>% rename(y = lpsa)
n <- nrow(train)
p <- ncol(train)

#--- BSS ----------

library(leaps)

fit_BSS <- regsubsets(y~.,train, nvmax=p)
summary_BSS <-summary(fit_BSS)
summary_BSS$outmat

fit_all <- regsubsets(y~.,train, nvmax=p, nbest=2^p, really.big = TRUE)
rss_all <- summary(fit_all)$rss
size_all <- apply(summary(fit_all)$which,1,sum)
rss_0 <- sum(summary(lm(y~1,train))$residuals^2)
#pdf("Figure_RSS_BSS.pdf")
plot(c(1,size_all)-1,c(rss_0,rss_all), xlab="Subset size", ylab="RSS", pch=19)
points(0:(p-1), c(rss_0,summary_BSS$rss), col=2, pch=19)
lines(0:(p-1), c(rss_0,summary_BSS$rss), col=2)
#dev.off()

#pdf("Figure_BIC_BSS.pdf")
plot(fit_BSS, scale="bic")
#dev.off()

#--- Forward Stepwise ----------

library(MASS)
fit_null <- lm(y ~1, data = train)
fit_full <- lm(y ~., data = train)
stepAIC(fit_null, 
        scope=list(upper=fit_full),
        direction = "forward", 
        k=0, trace = TRUE)

train_std <- data.frame(scale(train, center =TRUE, scale = sqrt(diag(var(train)*((n-1)/n)) ))[,])

fit_FS <- regsubsets(y~ .,train_std, intercept=FALSE, nvmax=p, method = "forward")
summary_FS <- summary(fit_FS)

RSQ <- summary_FS$rsq
beta_hat_mat <- matrix(0, nrow=p, ncol=p-1)
colnames(beta_hat_mat) <- names(train)[-9]
for (i in 1:(p-1) ){
beta_hat_mat[i+1, names(coef(fit_FS, i))]<- coef(fit_FS, i)
}
#pdf("Figure_R2_FS.pdf")
matplot(c(0,RSQ), beta_hat_mat, type="l", lwd=2, lty=1, xlab=expression(R^2), ylab=expression(hat(beta)))
abline(v=RSQ, lty=3)
#dev.off()

X_std = as.matrix(train_std[,-9])
y_std = train_std[,9]


#--- Lasso ----------

library(lars)
library(glmnet)


fit_lasso <- lars(x=X_std,y=y_std,type="lasso",intercept=FALSE,  normalize = FALSE)

#pdf("Figure_R2_Lasso.pdf")
matplot(fit_lasso$R2, fit_lasso$beta, type="l", lty=1, , ylab=expression(hat(beta)), xlab=expression(R^2), lwd=2)
abline(v=fit_lasso$R2[-1], lty=3)
axis(3, at=fit_lasso$R2, labels=c(round(fit_lasso$lambda,1),0))
#dev.off()




#--- Forward Stagewise ----------

forward_stagewise <- function(X,y,eps=0.01, itr = 100){
  
  n = nrow(X)
  p = ncol(X)
  r = y
  beta = rep(0,p)
  beta_mat <- matrix(0,ncol=p,nrow=itr)
  
  for (b in 1:itr){
  
  max_correlation = 0
  best_predictor = 0
  
  for (j in 1:p){
      current_correlation = max( abs(cor(r, X)) )
      if (current_correlation > max_correlation){
          best_predictor = which.max( abs(cor(r, X)) )
          max_correlation = current_correlation
      }
  }

  delta = eps * sign(crossprod(X[, best_predictor], r))
  beta[best_predictor] = beta[best_predictor] + delta
  
  beta_mat[b,] <- beta
  
  for (i in 1:n) r[i] = r[i] - delta * X[i, best_predictor]
  }
  return(beta_mat)
  
}

beta_hat_i <- forward_stagewise(X_std,y_std,eps=0.005, itr = 400)
L1_norm <- apply(beta_hat_i,1,function(x) sum(abs(x)))
#pdf("Figure_PATH_FSeps.pdf")
matplot(L1_norm, beta_hat_i, type="s", lty=1, xlab="L1 norm", ylab="Coefficients", lwd=2)
#dev.off()

fit_lasso <- glmnet(X_std,y_std,intercept=FALSE,standardize = FALSE)
#pdf("Figure_PATH_lasso.pdf")
plot(fit_lasso, lwd=2)
#dev.off()

X = as.matrix(train[,-9])
y = as.matrix(train[,9])

K <- 10
set.seed(123)
cv_fit <- cv.glmnet(X,y,nfolds = K)
#pdf("Figure_CV_lasso.pdf")
plot(cv_fit)
#dev.off()

cv_fit$lambda.1se
coef(cv_fit, s="lambda.1se")

cv_fit$lambda.min
coef(cv_fit, s="lambda.min")

set.seed(123)
fit_relax <- glmnet(X,y, relax = TRUE)

#pdf("Figure_relaxed_lasso.pdf")
plot(fit_relax, gamma = 0, lwd=2)
#dev.off()

set.seed(123)
cv_fit_relax <- cv.glmnet(X, y, relax = TRUE)
print(cv_fit_relax)
plot(cv_fit_relax, se.bands = FALSE)

cv_fit_relax0 <- cv.glmnet(X, y, gamma = 0, relax = TRUE)
#pdf("Figure_cv_relaxed0_lasso.pdf")
plot(cv_fit_relax0)
#dev.off()

rm(list=ls())

library(gglasso)

data(bardet)
group1 <- rep(1:20, each = 5)
fit_ls <- gglasso(x = bardet$x, y = bardet$y, group = group1, loss = "ls")
plot(fit_ls)
cvfit_ls <- cv.gglasso(x = bardet$x, y = bardet$y, group = group1, loss = "ls")
plot(cvfit_ls)
coef(cvfit_ls, s = "lambda.min")

#---------------------------------------
# HITTERS
#---------------------------------------

rm(list=ls())

library(tidymodels)
library(ISLR)

Hitters <- as_tibble(Hitters) %>%
  filter(!is.na(Salary))

set.seed(123)
Hitters_split <- initial_split(Hitters, strata = "Salary")

Hitters_train <- training(Hitters_split)
Hitters_test <- testing(Hitters_split)

library(leaps)
fit_BSS <- regsubsets(Salary~.,Hitters_train)
summary_BSS<-summary(fit_BSS)
head(summary_BSS$outmat, 10)

predict.regsubsets =function(object ,newdata ,id ,...){
  form=as.formula(object$call[[2]])
  mat=model.matrix(form, newdata)
  coefi =coef(object, id=id)
  xvars =names(coefi)
  mat[,xvars]%*%coefi
}

yhat_BSS_Cp = predict.regsubsets(fit_BSS, newdata=Hitters_test, id=which.min(summary_BSS$cp))
RMSE_BSS_Cp = sqrt(mean( (yhat_BSS_Cp - Hitters_test$Salary)^2 ))
RMSE_BSS_Cp

Hitters_fold <- vfold_cv(Hitters_train, v = 10)

elasticnet_recipe <- 
  recipe(formula = Salary ~ ., data = Hitters_train) %>% 
  step_novel(all_nominal_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors()) %>% 
  step_normalize(all_predictors())

elasticnet_spec <- 
  linear_reg(penalty = tune(), mixture = tune()) %>% 
  set_mode("regression") %>% 
  set_engine("glmnet") 

elasticnet_workflow <- workflow() %>% 
  add_recipe(elasticnet_recipe) %>% 
  add_model(elasticnet_spec)

tuning_grid <- grid_regular(parameters(penalty(range = c(-1, 2)), mixture(range = c(0, 1))), levels = c(50, 3))

tuning_grid %>% 
  ggplot(aes(x = mixture, y = penalty)) +
  geom_point() +
  scale_y_log10() +
  theme_bw()

tune_res <- tune_grid(
  elasticnet_workflow,
  resamples = Hitters_fold, 
  grid = tuning_grid
)

autoplot(tune_res) +
  theme_bw()

best_tuning <- select_best(tune_res, metric = "rmse")

elasticnet_final <- finalize_workflow(elasticnet_workflow, best_tuning)

elasticnet_final_fit <- fit(elasticnet_final, data = Hitters_train)

elasticnet_final_fit %>%
  extract_fit_engine() %>%
  plot(xvar = "lambda")

augment(elasticnet_final_fit, new_data = Hitters_test) %>%
  rmse(truth = Salary, estimate = .pred)
