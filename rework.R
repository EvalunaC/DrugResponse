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
#drug <- "Lapatinib"
#############################################
######### 0. Validation Function  #########
#############################################
methods_result <- function(drug_data, drug){
 trainFrame = drug_data
### RMSE and R2 function
eval_result<-function(preds){
  MAE <- mean(abs(trainFrame$Resp - preds))
  SSE <- sum((trainFrame$Resp - preds)^2)
  SST <- sum((trainFrame$Resp - mean(trainFrame$Resp))^2)
  SSM <- sum((preds-mean(trainFrame$Resp))^2)
  R_square <- 1 - SSE / SST
  R2 <- (cor(as.numeric(preds),trainFrame$Resp))^2
  R2adj<-1-((1-R_square)*(nrow(trainFrame)-1)/(nrow(trainFrame)-(ncol(trainFrame)-1)-1))
  RMSE = sqrt(SSE/nrow(trainFrame))
#  F_stat<-SSM/(ncol(trainFrame)-1)/(SSE/(nrow(trainFrame)-ncol(trainFrame)))
#  t_test<-t.test(trainFrame$Resp, preds)$p.value
#  ks_test<-ks.test(trainFrame$Resp, preds)$p.value
  results<-list(R2=R2,RMSE=RMSE,R_Square=R_square, Adjusted_R2=R2adj,MAE=MAE)
  return(results)

}
set.seed(1)

cat(paste(Sys.time(),"==========",drug,":\n"))
#############################################
######### 1. GR paper linear Ridge  #########
#############################################
cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Start...\n"))
model_GR <- linearRidge(Resp ~ ., data = trainFrame)
preds_GR<-predict(model_GR,trainFrame)
GR_result<-eval_result(preds_GR)
GR_result$method <- "GR paper linear Ridge"
GR_result$drug <- drug
cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Complete\n"))
#############################################
########## 2.  Random Forest
#############################################
cat(paste(Sys.time(),"==========","2. Random Forest Start...\n"))
model_rf <- train(Resp~.,data=trainFrame,
                         method = "rf",
                         ntree = 50,
                         tuneGrid = data.frame(mtry = 2),
                         nodesize = 5,
                         importance = TRUE,
                         metric = "RMSE",
                         trControl = trainControl(method = "oob", seed = c(1,1)),
                         allowParallel = FALSE)#save(model_rf,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_rf.RData")
preds_rf<-predict(model_rf,trainFrame)
rf_result<-eval_result(preds_rf)
rf_result$method <- "Random Forest"
rf_result$drug <- drug
cat(paste(Sys.time(),"==========","2. Random Forest Complete\n"))


#############################################
############## 3. Principle Component Regression
#############################################
cat(paste(Sys.time(),"==========","3. Principle Component Regression Start...\n"))
model_pcr<-train(Resp~.,data=trainFrame,method="pcr",importance=TRUE)
preds_pcr<-predict(model_pcr,trainFrame)
pcr_result<-eval_result(preds_pcr)
pcr_result$method <- "Principle Component Regression"
pcr_result$drug <- drug

cat(paste(Sys.time(),"==========","3. Principle Component Regression Complete\n"))
#top <- var_pcr$importance[with(var_pcr$importance,order(-Overall)),,drop = FALSE]
#var100_pcr <- top[1:100,,drop = FALSE]
#############################################
##############  4. Partial Least Square
#############################################
cat(paste(Sys.time(),"==========","4. Partial Least Square Start...\n"))
model_pls<-train(Resp~.,data=trainFrame,method="pls",importance=TRUE)
preds_pls<-predict(model_pls,trainFrame)
pls_result<-eval_result(preds_pls)
pls_result$method <- "Partial Least Square"
pls_result$drug <- drug
cat(paste(Sys.time(),"==========","4. Partial Least Square Complete\n"))
#############################################
############ 5. Ridge GLM penalty alpha=0
#############################################
cat(paste(Sys.time(),"==========","5. Ridge GLM penalty Start...\n"))
#model_ridgeglm<-train(Resp~.,data=trainFrame,method="glmnet",importance=TRUE,tuneGrid = expand.grid(alpha = 1)),trControl=trainControl("cv",number=10))
cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
(best_lam_0 <- cv_output_0$lambda.min) ###104.2796                             ########################### Evry time different. 49.54127, 47.28955, 114.4468
model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)  ######################Cannot run from my end.
#top <- var_ridgeglm[with(var_ridgeglm,order(-Overall)),,drop = FALSE]
#var100_ridgeglm <- top[1:100,,drop = FALSE]
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(trainFrame[,-1]))
RidgeGLM_result <-eval_result(preds_ridgeglm)
RidgeGLM_result$method <- "Ridge GLM"
RidgeGLM_result$drug <- drug
cat(paste(Sys.time(),"==========","5. Ridge GLM penalty Complete\n"))

#############################################
############ 6. Lasso GLM penalty alpha=1
#############################################
cat(paste(Sys.time(),"==========","6. Lasso GLM penalty Start...\n"))
cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
(best_lam_1 <- cv_output_1$lambda.min) ### 0.1444154
model_Lasso_1<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1, lambda=best_lam_1)
#top <- var_Lasso_1[with(var_Lasso_1,order(-Overall)),,drop = FALSE]
#var100_Lasso_1 <- top[1:100,,drop = FALSE]
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(trainFrame[,-1]))
Lasso_result_1<-eval_result(preds_Lasso_1)
Lasso_result_1$method <- "Lasso GLM"
Lasso_result_1$drug <- drug
cat(paste(Sys.time(),"==========","6. Lasso GLM penalty Complete\n"))

#############################################
############ 7. K-nearest neighbors (KNN) algorithm
#############################################
cat(paste(Sys.time(),"==========","7. K-nearest neighbors (KNN) algorithm Start...\n"))
model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
preds_knn<-predict(model_KNN,trainFrame)
knn_result<-eval_result(preds_knn)
knn_result$method <- "K-nearest neighbors (KNN) algorithm"
knn_result$drug <- drug
cat(paste(Sys.time(),"==========","7. K-nearest neighbors (KNN) algorithm Complete\n"))
#############################################
############ 8. Support vector machine regression
#############################################
cat(paste(Sys.time(),"==========","8. Support vector machine regression Start...\n"))
model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2',importance=TRUE)
preds_svm<-predict(model_svm,trainFrame)
svm_result<-eval_result(preds_svm)

svm_result$method <- "Support vector machine regression"
svm_result$drug <- drug

cat(paste(Sys.time(),"==========","8. Support vector machine regression Complete\n"))

#############################################
############ 9. Treebag (bootstrap aggregating) algorithm
#############################################
cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Start...\n"))
#model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag')
#preds_treebag<-predict(model_treebag,trainFrame)
#treebag_result<-eval_result
#(preds_treebag)
#treebag_result$method <- "Treebag (bootstrap aggregating) algorithm"
#treebag_result$drug <- drug
cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Complete\n"))

#############################################
############ 10. Elastic Net
#############################################
cat(paste(Sys.time(),"==========","10. Elastic Net Regression Start...\n"))
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))

preds_EN<-predict(model_EN,trainFrame)
en_result<-eval_result(preds_EN)
en_result$method <- "Elastic Net"
en_result$drug <- drug
cat(paste(Sys.time(),"==========","10. Elastic Net Regression Complete\n"))

#########################################################################################
#########################            combine results               ######################
#########################################################################################
l <- list(GR_result,rf_result,pcr_result,pls_result,RidgeGLM_result,Lasso_result_1,knn_result,svm_result,en_result)
#l <- list(GR_result,rf_result,pcr_result,pls_result,RidgeGLM_result,Lasso_result_1,knn_result,svm_result,treebag_result,en_result)
Result_final <- ldply(l, data.frame)
df <- Result_final[,c("method","drug","RMSE","R_Square","MAE","R2")]
  name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output_04092021/",drug,".csv", sep="",collapse = NULL)
  write.csv(df,name, quote=FALSE)
}
