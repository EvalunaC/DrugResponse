source("Creating_drug_data.R")
memory.limit(9999999999)

library(caret)
library(pRRophetic)
library(ridge)
library(MASS)
library(dplyr)
library(caret)
library(glmnet)
library(randomForest)
library(pls)
library(plyr)
library(deepnet)

library(ranger)
library(kknn)
library(rfinterval)

set.seed(1)


cat(paste(Sys.time(),"==========",drug,":\n"))

cat(paste(Sys.time(),"==========","GR paper linear Ridge Start...\n"))
model_GR <- linearRidge(Resp ~ ., data = trainFrame)

cat(paste(Sys.time(),"==========","Random Forest Start...\n"))
model_rf <- train(Resp~.,data=trainFrame,
                         method = "rf",
                         ntree = 50,
                         tuneGrid = data.frame(mtry = 2),
                         nodesize = 5,
                         importance = TRUE,
                         metric = "RMSE",
                         trControl = trainControl(method = "oob", seed = c(1,1)),
                         allowParallel = FALSE)

cat(paste(Sys.time(),"==========","Random Forest Ranger Start...\n"))
#model_ranger<-train(Resp~.,data=trainFrame,num.trees = 50,method="ranger",trControl = trainControl(method = "oob", seed = list(c(1,1,1,1,1,1),c(1,1,1,1,1,1))))
model_ranger<-ranger(Resp~.,data=trainFrame,num.trees = 50,splitrule="variance",seed =1)


cat(paste(Sys.time(),"==========","Principle Component Regression Start...\n"))
model_pcr<-train(Resp~.,data=trainFrame,method="pcr")#,importance=TRUE)

cat(paste(Sys.time(),"==========","Partial Least Square Start...\n"))
model_pls<-train(Resp~.,data=trainFrame,method="pls")#,importance=TRUE)

cat(paste(Sys.time(),"==========","K-nearest neighbors (KNN) algorithm Start...\n"))
model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
cat(paste(Sys.time(),"==========","Weighted K-nearest neighbors (KNN) algorithm Start...\n"))

model_KKNN<-train(Resp~.,data=trainFrame,method="kknn",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","Support vector machine regression Start...\n"))
model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2')#,importance=TRUE)

cat(paste(Sys.time(),"==========","Treebag (bootstrap aggregating) algorithm Start...\n"))
model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag')
#preds_treebag<-predict(model_treebag,trainFrame)
#treebag_result<-eval_result
#(preds_treebag)
#treebag_result$method <- "Treebag (bootstrap aggregating) algorithm"
#treebag_result$drug <- drug
#cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Complete\n"))

cat(paste(Sys.time(),"==========","Elastic Net Regression Start...\n"))
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","Ridge GLM penalty Start...\n"))
cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
(best_lam_0 <- cv_output_0$lambda.min) ###104.2796    Evry time different. 49.54127, 47.28955, 114.4468
model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)

cat(paste(Sys.time(),"==========","Lasso GLM penalty Start...\n"))
cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
(best_lam_1 <- cv_output_1$lambda.min) ### 0.1444154
model_Lasso_1<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1, lambda=best_lam_1)

#### Regression-Enhanced Random Forests, "Predictive Inference for Random Forests" RF+Lasso
cat(paste(Sys.time(),"==========","Predictive Inference for Random Forests (oob) Start...\n"))

model_rfinv1<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
     method = "oob", alpha = 0.05,
     symmetry = TRUE)

cat(paste(Sys.time(),"==========","Predictive Inference for Random Forests (split-conformal) Start...\n"))
model_rfinv2<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
    method = "split-conformal", alpha = 0.05,
    seed = 1)
cat(paste(Sys.time(),"==========","Predictive Inference for Random Forests (quantreg) Start...\n"))
model_rfinv3<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
    method = "quantreg", alpha = 0.05)

#model_rfinv3$testPred
