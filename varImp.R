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

### RMSE and R2 function
eval_result<-function(trainFrame,preds){
  MAE <- mean(abs(trainFrame$Resp - preds))
  SSE <- sum((trainFrame$Resp - preds)^2)
  SST <- sum((trainFrame$Resp - mean(trainFrame$Resp))^2)
  SSM <- sum((preds-mean(trainFrame$Resp))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(trainFrame))
  results<-list(RMSE=RMSE,R_Square=R_square, MAE=MAE)
  return(results)
}



methods_result <- function(drug_data, drug){
 trainFrame = drug_data
cat(paste(Sys.time(),"==========",drug,":\n"))
#############################################
######### 1. GR paper linear Ridge  #########
#############################################
cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Start...\n"))
model_GR <- linearRidge(Resp ~ ., data = trainFrame)
preds_GR<-predict(model_GR,trainFrame)
GR_result<-eval_result(trainFrame,preds_GR)
GR_result$method <- "GR paper linear Ridge"
GR_result$drug <- drug
cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Complete\n"))
#############################################
########## 2.  Random Forest
#############################################
cat(paste(Sys.time(),"==========","2. Random Forest Start...\n"))
model_rf<-train(Resp~.,data=trainFrame,method="rf",importance=TRUE)
#save(model_rf,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_rf.RData")
best <- as.numeric(row.names(model_rf$bestTune))
rf_result <- list(method = "Random Forest",
                   drug = drug,
                   RMSE = model_rf$result$RMSE[best],
                   R_Square = model_rf$result$Rsquared[best],
                   MAE = model_rf$results$MAE[best]
                  )
cat(paste(Sys.time(),"==========","2. Random Forest Complete\n"))
var_rf <- varImp(model_rf ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_rf$importance[with(var_rf$importance,order(-Overall)),,drop = FALSE]
#var100_rf <- top[1:100,,drop = FALSE]
#############################################
############## 3. Principle Component Regression
#############################################
cat(paste(Sys.time(),"==========","3. Principle Component Regression Start...\n"))
#model_pcr<-pcr(Resp~.,data=trainFrame,ncomp=3, validation = "CV", jackknife = TRUE)
model_pcr<-train(Resp~.,data=trainFrame,method="pcr",importance=TRUE)
#save(model_pcr,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_pcr.RData")
best <- as.numeric(row.names(model_pcr$bestTune))
pcr_result <- list(method = "Principle Component Regression",
                   drug = drug,
                   RMSE = model_pcr$result$RMSE[best],
                   R_Square = model_pcr$result$Rsquared[best],
                   MAE = model_pcr$results$MAE[best]
                  )
cat(paste(Sys.time(),"==========","3. Principle Component Regression Complete\n"))
var_pcr <- varImp(model_pcr ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_pcr$importance[with(var_pcr$importance,order(-Overall)),,drop = FALSE]
#var100_pcr <- top[1:100,,drop = FALSE]
#############################################
##############  4. Partial Least Square
#############################################


cat(paste(Sys.time(),"==========","4. Partial Least Square Start...\n"))
model_pls<-train(Resp~.,data=trainFrame,method="pls",importance=TRUE)
#save(model_pls,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_pls.RData")
best <- as.numeric(row.names(model_pls$bestTune))
pls_result <- list(method = "Partial Least Square",
                   drug = drug,
                   RMSE = model_pls$result$RMSE[best],
                   R_Square = model_pls$result$Rsquared[best],
                   MAE = model_pls$results$MAE[best]
                  )
var_pls <- varImp(model_pls ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_pls$importance[with(var_pls$importance,order(-Overall)),,drop = FALSE]
#var100_pls <- top[1:100,,drop = FALSE]
cat(paste(Sys.time(),"==========","4. Partial Least Square Complete\n"))
#############################################
############ 5. Ridge GLM penalty alpha=0
#############################################
cat(paste(Sys.time(),"==========","5. Ridge GLM penalty Start...\n"))
#model_ridgeglm<-train(Resp~.,data=trainFrame,method="glmnet",importance=TRUE,tuneGrid = expand.grid(alpha = 1)),trControl=trainControl("cv",number=10))
cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
(best_lam_0 <- cv_output_0$lambda.min) ###104.2796                             ########################### Evry time different. 49.54127, 47.28955, 114.4468
model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)  ######################Cannot run from my end.
var_ridgeglm <- varImp(model_ridgeglm ,useModel = TRUE, nonpara = TRUE, scale = TRUE,lambda=best_lam_0)
#top <- var_ridgeglm[with(var_ridgeglm,order(-Overall)),,drop = FALSE]
#var100_ridgeglm <- top[1:100,,drop = FALSE]
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(trainFrame[,-1]))
RidgeGLM_result <-eval_result(trainFrame,preds_ridgeglm)
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
var_Lasso_1 <- varImp(model_Lasso_1 ,useModel = TRUE, nonpara = TRUE, scale = TRUE,lambda=best_lam_0)
#top <- var_Lasso_1[with(var_Lasso_1,order(-Overall)),,drop = FALSE]
#var100_Lasso_1 <- top[1:100,,drop = FALSE]
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(trainFrame[,-1]))
Lasso_result_1<-eval_result(trainFrame,preds_Lasso_1)
Lasso_result_1$method <- "Lasso GLM"
Lasso_result_1$drug <- drug
cat(paste(Sys.time(),"==========","6. Lasso GLM penalty Complete\n"))

#############################################
############ 7. K-nearest neighbors (KNN) algorithm
#############################################
cat(paste(Sys.time(),"==========","7. K-nearest neighbors (KNN) algorithm Start...\n"))
model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
#save(model_KNN,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_KNN.RData")
best <- as.numeric(row.names(model_KNN$bestTune))
KNN_result <- list(method = "K-nearest neighbors (KNN) algorithm",
                   drug = drug,
                   RMSE = model_KNN$result$RMSE[best],
                   R_Square = model_KNN$result$Rsquared[best],
                   MAE = model_KNN$results$MAE[best]
                  )
var_KNN <- varImp(model_KNN ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_KNN$importance[with(var_KNN$importance,order(-Overall)),,drop = FALSE]
#var100_KNN <- top[1:100,,drop = FALSE]
cat(paste(Sys.time(),"==========","7. K-nearest neighbors (KNN) algorithm Complete\n"))
#############################################
############ 8. Support vector machine regression
#############################################
cat(paste(Sys.time(),"==========","8. Support vector machine regression Start...\n"))
model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2',importance=TRUE)
#save(model_svm,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_svm.RData")
#preds_svm<-predict(model_svm,trainFrame)
#svm_result<-eval_result(trainFrame,preds_svm)
best <- as.numeric(row.names(model_svm$bestTune))
svm_result <- list(method = "Support vector machine regression",
                   drug = drug,
                   RMSE = model_svm$result$RMSE[best],
                   R_Square = model_svm$result$Rsquared[best],
                   MAE = model_svm$results$MAE[best]
                  )
var_svm <- varImp(model_svm ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_svm$importance[with(var_svm$importance,order(-Overall)),,drop = FALSE]
#var100_svm <- top[1:10,,drop = FALSE]

cat(paste(Sys.time(),"==========","8. Support vector machine regression Complete\n"))

#############################################
############ 9. Treebag (bootstrap aggregating) algorithm
#############################################
cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Start...\n"))
model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag')
#save(model_treebag,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_treebag.RData")
best <- as.numeric(row.names(model_treebag$bestTune))
treebag_result <- list(method = "Treebag (bootstrap aggregating) algorithm",
                   drug = drug,
                   RMSE = model_treebag$result$RMSE[best],
                   R_Square = model_treebag$result$Rsquared[best],
                   MAE = model_treebag$results$MAE[best]
                  )
var_treebag <- varImp(model_treebag ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_treebag$importance[with(var_treebag$importance,order(-Overall)),,drop = FALSE]
#var100_treebag <- top[1:10,,drop = FALSE]
cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Complete\n"))

#############################################
############ 10. Elastic Net
#############################################
cat(paste(Sys.time(),"==========","10. Elastic Net Regression Start...\n"))
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))
#save(model_EN,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_EN.RData")
best <- as.numeric(row.names(model_EN$bestTune))
EN_result <- list(method = "Elastic Net",
                   drug = drug,
                   RMSE = model_EN$result$RMSE[best],
                   R_Square = model_EN$result$Rsquared[best],
                   MAE = model_EN$results$MAE[best]
                  )

var_EN <- varImp(model_EN ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
#top <- var_EN$importance[with(var_EN$importance,order(-Overall)),,drop = FALSE]
#var10_EN <- top[1:10,,drop = FALSE]
cat(paste(Sys.time(),"==========","10. Elastic Net Regression Complete\n"))

#########################################################################################
#########################            combine results               ######################
#########################################################################################
l <- list(GR_result, RidgeGLM_result,rf_result, pcr_result, pls_result,Lasso_result_1,svm_result,treebag_result,EN_result,KNN_result)

Result_final <- ldply(l, data.frame)
df <- Result_final[,c("method","drug","RMSE","R_Square","MAE")]
  name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output_04092021/",drug,".csv", sep="",collapse = NULL)
  write.csv(df,name, quote=FALSE)

#########################################################################################
######################### variable importance evaluation functions ######################
#########################################################################################
top <- var_rf$importance[with(var_rf$importance,order(-Overall)),,drop = FALSE]
var10_rf <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),rf = seq(1,10))

top <- var_pcr$importance[with(var_pcr$importance,order(-Overall)),,drop = FALSE]
var10_pcr <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),pcr = seq(1,10))

top <- var_pls$importance[with(var_pls$importance,order(-Overall)),,drop = FALSE]
var10_pls <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),pls = seq(1,10))

top <- var_ridgeglm[with(var_ridgeglm,order(-Overall)),,drop = FALSE]
var10_ridgeglm <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),ridgeglm = seq(1,10))

top <- var_Lasso_1[with(var_Lasso_1,order(-Overall)),,drop = FALSE]
var10_Lasso_1 <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),Lasso_1 = seq(1,10))

top <- var_KNN$importance[with(var_KNN$importance,order(-Overall)),,drop = FALSE]
var10_KNN <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),KNN = seq(1,10))

top <- var_svm$importance[with(var_svm$importance,order(-Overall)),,drop = FALSE]
var10_svm <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),svm = seq(1,10))

top <- var_treebag$importance[with(var_treebag$importance,order(-Overall)),,drop = FALSE]
var10_treebag <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),treebag = seq(1,10))

top <- var_EN$importance[with(var_EN$importance,order(-Overall)),,drop = FALSE]
var10_EN <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),EN = seq(1,10))

data1 <- merge(var10_rf,var10_pcr,by="gene",all = T,suffixes = c("",""))
data2 <- merge(var10_pls,var10_ridgeglm,by="gene",all = T,suffixes = c("",""))
data3 <- merge(var10_Lasso_1,var10_KNN,by="gene",all = T,suffixes = c("",""))
data4 <- merge(var10_svm,var10_treebag,by="gene",all = T,suffixes = c("",""))

data5 <- merge(data1,data2,by="gene",suffixes = c("",""),all = T)
data6 <- merge(data3,data4,by="gene",suffixes = c("",""),all = T)
data7 <- merge(data5,data6,by="gene",suffixes = c("",""),all = T)
data <- merge(data7,var10_EN,by="gene",suffixes = c("",""),all = T)

data$overlap <- rowSums(!is.na(data))-1
imp_rank <- data[order(-data$overlap),]

name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output_varImp_04092021/",drug,".csv", sep="",collapse = NULL)
write.csv(imp_rank,name, quote=FALSE)
}
