## This script runs classifiers to learn patterns associated with matched / unmatched records.
##
## Author: Yoni Friedman, Palladium
## Last Edited: April 11, 2021

# Load libraries and data ----------------------
# Set Working Directory
setwd("~/Kenya/Deduplication")

# Load Libraries
library(dplyr)
library(ggplot2)
library(caret)
library(randomForest)
library(xgboost)
library(reshape2)
library(PRROC)
library(neuralnet)
library(lime)

# Load Data
dat <- readRDS('./matches_prep.rds')

# Look for collinearity -----------------------

# generate correlation matrix
cormat <- cor(dat)

# Identify any variables with perfect correlation
cormat_long <- melt(cormat) # Get into long format
head(cormat_long)
cormat_perf <- cormat_long %>%
  filter(value > .99) %>%
  filter(Var1 != Var2)

# Drop variables collinear with other variables
cols_to_drop <- c("fnSXLV", "mnSXLV", "lnSXLV", "YOBLV", "MOBLV", "DOBLV")
dat <- dat[, which(!names(dat) %in% cols_to_drop)]

# Split into train / eval / test ---------------------

# Will set up a 60/20/20 train/eval/test split

# Set seed for reproducability
set.seed(2231)

# Split 80/20 between train and test
split <- createDataPartition(y = dat$target,p = 0.8,list = FALSE)
train_all <- dat[split, ] 
test <- dat[-split, ] # hold out

# Split train_all into train and eval by 75/25, which is 60/20 of the overall data
split_val <- createDataPartition(y = train_all$target,p = 0.75,list = FALSE)
train <- train_all[split_val, ]
val <- train_all[-split_val, ]

# Train logistic regression ----------------------

log_model <- glm(target ~ ., data = train,
                 family = binomial(link="logit"))

summary(log_model)

log_predict <- predict(log_model,newdata = val,type = "response")
fg <- log_predict[val$target == 1]
bg <- log_predict[val$target == 0]

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr)

# Train Random Forest ---------------------------
grid <- expand.grid(mtry = c(4,6,8), nodesize = c(5, 1))

for (i in 1:nrow(grid)){
  set.seed(2231)
  
  rf <- randomForest(
    as.factor(target) ~ .,
    data = train,
    mtry = grid[i, 1],
    nodesize = grid[i, 2],
    importance = TRUE
  )
  
  pred_val = predict(rf, newdata=val, type = "prob")
  fg <- pred_val[val$target == 1, 2]
  bg <- pred_val[val$target == 0, 2]
  prc <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  grid$val_pr_auc[i] <- prc$auc.integral
  
}

# Print output
grid %>% arrange(desc(val_pr_auc))

# Train best model
set.seed(2231)

rf <- randomForest(
  as.factor(target) ~ .,
  data = train,
  mtry = 6,
  nodesize = 1,
  importance = TRUE
)

# Plot Feature Importance
varImpPlot(rf,type=2)

# Train XGBoost ---------------------------------
grid <- expand.grid(eta = c(0.1, 0.3),
                    max_depth = c(4, 6, 8, 10),
                    cs = c(.3, .5, .7, .9))

for (i in 1:nrow(grid)){
  
  set.seed(2231)
  xgb <- xgboost::xgboost(data = data.matrix(train[,-1]), 
                          label = train$target, 
                          eta = grid[i, 1],
                          max_depth = grid[i, 2], 
                          nround=25, 
                          # subsample = 0.5,
                          colsample_bytree = grid[i, 3],
                          seed = 2231,
                          objective = "binary:logistic",
                          metric = 'auc',
                          # stratified = TRUE,
                          verbose = 0
  )
  
  val_predict <- predict(xgb,newdata = data.matrix(val[, -1]))
  fg <- val_predict[val$target == 1]
  bg <- val_predict[val$target == 0]
  prc <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  grid$val_pr_auc[i] <- prc$auc.integral
  
}

grid %>% arrange(desc(val_pr_auc))

# Train best model
set.seed(2231)
xgb <- xgboost::xgboost(data = data.matrix(train[,-1]), 
                        label = train$target, 
                        eta = .3,
                        max_depth = 10, 
                        nround=25, 
                        # subsample = 0.5,
                        colsample_bytree = .9,
                        seed = 2231,
                        objective = "binary:logistic",
                        metric = 'auc',
                        # stratified = TRUE,
                        verbose = 0
)

# Plot feature importance
mat <- xgb.importance(feature_names = colnames(train[,-1]),model = xgb)
xgb.plot.importance(importance_matrix = mat[1:20]) #first 20 variables

# Train Neural Net -------------------------------

# Set seed for reproducability
set.seed(2231)

# Set grid for three layer neural network
grid <- expand.grid(layer_1 = c(20, 10),
                    layer_2 = c(10),
                    layer_3 = c(10, 5))

for (i in 1:nrow(grid)){
  
  set.seed(2231)
  NN <- neuralnet(target ~ .,
                  data = train,
                  hidden = c(grid[i, 1], grid[i, 2], grid[i, 3]),
                  linear.output = FALSE
  )
  
  val_predict = compute(NN, val[,c(2:ncol(dat))])$net.result
  fg <- val_predict[val$target == 1]
  bg <- val_predict[val$target == 0]
  prc <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  grid$val_pr_auc[i] <- prc$auc.integral
  
}

grid %>% arrange(desc(val_pr_auc))
