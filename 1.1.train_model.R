set.seed(1000)
memory.limit(9999999999)

#setwd("/extraspace/ychen42/Drug_Response/Data/")
#load("brcaTrainFrame.RData")
#load("C:/Users/qiangli/Desktop/Drug Response Project/DR_Pred/trainFrame.RData")
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
set.seed(1)

#train_model <- function(drug,trainFrame,out_model_file){
#
#method_list <- c("GR paper linear Ridge",
#            "Random Forest",
#            "Principle Component Regression",
#            "Partial Least Square",
#            "KNN",
#            "SVM",
##           "Tree Bag",
#            "Elastic Net",
#            "Ridge GLM",
#            "Lasso GLM"
#            )
#
cat(paste(Sys.time(),"==========",drug,":\n"))

cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Start...\n"))
model_GR <- linearRidge(Resp ~ ., data = trainFrame)

cat(paste(Sys.time(),"==========","2. Random Forest Start...\n"))
model_rf <- train(Resp~.,data=trainFrame,
                         method = "rf",
                         ntree = 50,
                         tuneGrid = data.frame(mtry = 2),
                         nodesize = 5,
                         importance = TRUE,
                         metric = "RMSE",
                         trControl = trainControl(method = "oob", seed = c(1,1)),
                         allowParallel = FALSE)


model_ranger<-train(Resp~.,data=trainFrame,num.trees = 50,method="ranger",trControl = trainControl(method = "oob", seed = c(1,1)))

cat(paste(Sys.time(),"==========","3. Principle Component Regression Start...\n"))
model_pcr<-train(Resp~.,data=trainFrame,method="pcr",importance=TRUE)

cat(paste(Sys.time(),"==========","4. Partial Least Square Start...\n"))
model_pls<-train(Resp~.,data=trainFrame,method="pls",importance=TRUE)

cat(paste(Sys.time(),"==========","5. K-nearest neighbors (KNN) algorithm Start...\n"))
model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))

model_KKNN<-train(Resp~.,data=trainFrame,method="kknn",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","6. Support vector machine regression Start...\n"))
model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2',importance=TRUE)

#############################################
############ 9. Treebag (bootstrap aggregating) algorithm
#############################################
#cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Start...\n"))
#model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag')
#preds_treebag<-predict(model_treebag,trainFrame)
#treebag_result<-eval_result
#(preds_treebag)
#treebag_result$method <- "Treebag (bootstrap aggregating) algorithm"
#treebag_result$drug <- drug
#cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Complete\n"))

cat(paste(Sys.time(),"==========","7. Elastic Net Regression Start...\n"))
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","8. Ridge GLM penalty Start...\n"))
cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
(best_lam_0 <- cv_output_0$lambda.min) ###104.2796    Evry time different. 49.54127, 47.28955, 114.4468
model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)

cat(paste(Sys.time(),"==========","9. Lasso GLM penalty Start...\n"))
cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
(best_lam_1 <- cv_output_1$lambda.min) ### 0.1444154
model_Lasso_1<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1, lambda=best_lam_1)
