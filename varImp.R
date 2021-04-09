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
  R2adj<-1-((1-R_square)*(nrow(trainFrame)-1)/(nrow(trainFrame)-(ncol(trainFrame)-1)-1))
  RMSE = sqrt(SSE/nrow(trainFrame))
  F_stat<-SSM/(ncol(trainFrame)-1)/(SSE/(nrow(trainFrame)-ncol(trainFrame)))
  t_test<-t.test(trainFrame$Resp, preds)$p.value
  ks_test<-ks.test(trainFrame$Resp, preds)$p.value
  results<-list(RMSE=RMSE,R_Square=R_square, Adjusted_R2=R2adj,MAE=MAE, F_stat=F_stat,t_test=t_test,ks_test=ks_test)
  return(results)
}

##### AIC and BIC function for glmnet result

AIC_BIC_glmnet<-function(model){
  tLL <- model$nulldev - deviance(model)
  k <- model$df
  n <- model$nobs
  AIC <- -tLL+2*k
  AICc <- - tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
  BIC<-log(n)*k - tLL
  return(list(AIC=AIC,AICc=AICc,BIC=BIC))
}

#### AIC and BIC function for tain result
AIC_BIC_train<-function(model){
  tLL <- model$finalModel$nulldev - model$finalModel$nulldev * (1 -  model$finalModel$dev.ratio)
  k <- model$finalModel$df
  n <- model$finalModel$nobs
  AIC <- -tLL+2*k
  AICc <- - tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1)
  BIC<-log(n)*k - tLL
  return(list(AIC=AIC,AICc=AICc,BIC=BIC))
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

GR_result<-eval_result(trainFrame,preds_GR) ## RMSE=0.639 R2=0.729 R2_adj=1.01 MAE=0.499 AIC=

GR_result$method <- "GR paper linear Ridge"
GR_result$drug <- drug

varImp(model_GR)

cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Complete\n"))

var_GR <- varImp(model_GR)

#############################################
########## 2.  Random Forest
#############################################
cat(paste(Sys.time(),"==========","2. Random Forest Start...\n"))

system.time(model_rf<-train(Resp~.,data=trainFrame,method="rf",importance=TRUE))

save(model_rf,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_rf.RData")

best <- as.numeric(row.names(model_rf$bestTune))
rf_result <- list(method = "Random Forest",
                   drug = drug,
                   RMSE = model_rf$result$RMSE[best],
                   R_Square = model_rf$result$Rsquared[best],
                   MAE = model_rf$results$MAE[best]
                  )

cat(paste(Sys.time(),"==========","2. Random Forest Complete\n"))
var_rf <- varImp(model_rf ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
top <- var_rf$importance[with(var_rf$importance,order(-Overall)),,drop = FALSE]
var100_rf <- top[1:100,,drop = FALSE]

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
top <- var_pcr$importance[with(var_pcr$importance,order(-Overall)),,drop = FALSE]
var100_pcr <- top[1:100,,drop = FALSE]

#############################################
##############  4. Partial Least Square
#############################################
cat(paste(Sys.time(),"==========","4. Partial Least Square Start...\n"))
model_pls<-train(Resp~.,data=trainFrame,method="pls",importance=TRUE)
save(model_pls,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_pls.RData")


best <- as.numeric(row.names(model_pls$bestTune))
pls_result <- list(method = "Partial Least Square",
                   drug = drug,
                   RMSE = model_pls$result$RMSE[best],
                   R_Square = model_pls$result$Rsquared[best],
                   MAE = model_pls$results$MAE[best]
                  )


var_pls <- varImp(model_pls ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
top <- var_pls$importance[with(var_pls$importance,order(-Overall)),,drop = FALSE]
var100_pls <- top[1:100,,drop = FALSE]


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
top <- var_ridgeglm[with(var_ridgeglm,order(-Overall)),,drop = FALSE]
var100_ridgeglm <- top[1:100,,drop = FALSE]


preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(trainFrame[,-1]))

RidgeGLM_result <-eval_result(trainFrame,preds_ridgeglm)
RidgeGLM_result$AIC<-AIC_BIC_glmnet(model_ridgeglm)$AIC
RidgeGLM_result$AICc<-AIC_BIC_glmnet(model_ridgeglm)$AICc
RidgeGLM_result$BIC<-AIC_BIC_glmnet(model_ridgeglm)$BIC

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
top <- var_Lasso_1[with(var_Lasso_1,order(-Overall)),,drop = FALSE]
var100_Lasso_1 <- top[1:100,,drop = FALSE]



preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(trainFrame[,-1]))

Lasso_result_1<-eval_result(trainFrame,preds_Lasso_1)
Lasso_result_1$AIC<-AIC_BIC_glmnet(model_Lasso_1)$AIC
Lasso_result_1$AICc<-AIC_BIC_glmnet(model_Lasso_1)$AICc
Lasso_result_1$BIC<-AIC_BIC_glmnet(model_Lasso_1)$BIC

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
top <- var_KNN$importance[with(var_KNN$importance,order(-Overall)),,drop = FALSE]
var100_KNN <- top[1:100,,drop = FALSE]

cat(paste(Sys.time(),"==========","7. K-nearest neighbors (KNN) algorithm Complete\n"))

#############################################
############ 8. Support vector machine regression
#############################################

cat(paste(Sys.time(),"==========","8. Support vector machine regression Start...\n"))

model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2',importance=TRUE)
save(model_svm,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_svm.RData")

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
top <- var_svm$importance[with(var_svm$importance,order(-Overall)),,drop = FALSE]
var100_svm <- top[1:10,,drop = FALSE]



#svm_result$method <- "Support vector machine regression"
#svm_result$drug <- drug
#svm_result$RMSE <- model_svm$result$RMSE[3]
#svm_result$Rsquared <- model_svm$result$Rsquared[3]
#svm_result$MAE <- model_svm$results$MAE[3]
cat(paste(Sys.time(),"==========","8. Support vector machine regression Complete\n"))

#svm_result$RMSESD <- model_svm$results[2,5]
#svm_result$RsquaredSD <- model_svm$results[2,6]
#svm_result$MAESD <- model_svm$results[2,7]

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
top <- var_treebag$importance[with(var_treebag$importance,order(-Overall)),,drop = FALSE]
var100_treebag <- top[1:10,,drop = FALSE]



cat(paste(Sys.time(),"==========","9. Treebag (bootstrap aggregating) algorithm Complete\n"))

#############################################
############ 10. Elastic Net
#############################################
cat(paste(Sys.time(),"==========","10. Elastic Net Regression Start...\n"))
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))
save(model_EN,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_EN.RData")


best <- as.numeric(row.names(model_EN$bestTune))
EN_result <- list(method = "Elastic Net",
                   drug = drug,
                   RMSE = model_EN$result$RMSE[best],
                   R_Square = model_EN$result$Rsquared[best],
                   MAE = model_EN$results$MAE[best]
                  )

var_EN <- varImp(model_EN ,useModel = TRUE, nonpara = TRUE, scale = TRUE)
top <- var_EN$importance[with(var_EN$importance,order(-Overall)),,drop = FALSE]
var10_EN <- top[1:10,,drop = FALSE]



cat(paste(Sys.time(),"==========","10. Elastic Net Regression Complete\n"))







####################################################
############# 11. Stacked AutoEncoder Deep Neural Network
####################################################

#model_dnn<-train(Resp~.,data=trainFrame,method = 'dnn')
#save(model_dnn,file="/extraspace/ychen42/Drug_Response/yiqings_work/Lapatinib_models/model_dnn.RData")


#dnn_result <- list(method = "Stacked AutoEncoder Deep Neural Network",
#                   drug = drug,
#                   RMSE = model_dnn$result$RMSE,
#                   R_Square = model_dnn$result$Rsquared，
#                   MAE = model_dnn$result$MAE
#                  )


#l <- list(GR_result, RidgeGLM_result,RF_result, PCR_result, PLSR_result,Lasso_result_1,SVM_result,treebag_result,EN_result,KNN_result)
#l <- list(GR_result, RidgeGLM_result,RF_result, PCR_result, PLSR_result,Lasso_result_1,svm_result,treebag_result,EN_result,KNN_result)
#Result_final <- ldply(l, data.frame)
#df <- Result_final[,c("method","drug","RMSE","R_Square","Adjusted_R2","MAE","F_stat","t_test","ks_test","AIC","AICc","BIC")]
#  name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output/",drug,".csv", sep="",collapse = NULL)
#  write.csv(df,name, quote=FALSE)
#write.csv(df, "/extraspace/ychen42/Drug_Response/yiqings_work/Output/192Drug_MethodsResult.csv",
#          col.names = FALSE,
#          row.names = TRUE,
#          quote = FALSE)
#  return(df)





}

#drug <- "Vinblastine"

#name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output/",drug,".csv", sep="",collapse = NULL)
#write.csv(df,name, quote=FALSE)






#################### the following methods generate errors or taking forever

########### Robust Statistics  ## can not run
#robfit <- glmrob(Resp~.,data=trainFrame,family=gaussian,method="Mqle")
#plot(robfit)
#mean(robfit$residuals^2)
#summary(robfit)


##########################################
############## Quantile Regression Neural Network ## takes too long time
##########################################
library(qrnn)
model_qrnn<-train(Resp~.,data=trainFrame,method = 'qrnn',trControl=trainControl("cv",number=10))
preds_qrnn<-predict(model_qrnn,trainFrame)

model_qrnn <- qrnn.fit(x=as.matrix(trainFrame[,2:13543]),y=as.matrix(trainFrame$Resp),data=trainFrame,method="adam",n.hidden=2, n.trials=1, tau=0.2678857,n.errors.max=1000,iter.max=50,bag=TRUE) ## Error

############## ############## ##############
############## ???? Robust Linear Model
############## ############## ##############
model_rlm<-train(Resp~.,data=trainFrame,method = 'rlm'),trControl=trainControl("cv",number=10)) ## Error
preds_rlm<-predict(model_rlm,testFrame)
#preds_rlm<- preds_rlm^(1/0.8944584) # no boxcox for now.
#preds_rlm <- preds_rlm - 2.944905

model_rlm<-rlm(Resp~.,data=trainFrame,psi = psi.bisquare) ##Error: x is singular

#################################################################
#############  Quantile Regression - conquer method
##############################################################################
library(quantreg)
 model_rq <- rq(Resp~.,data=trainFrame, method = "conquer") the number of conlumns of x can not exceed the number of rows x

library(conquer)
model_conquer <- conquer(X=as.matrix(trainFrame[,-1]),Y=trainFrame$Resp,checkSing = TRUE) #same with above

####################################################
############# Bayesian Regularized Neural Networks
####################################################
library(brnn)
library(parallel)
model_brnn<-brnn(Resp~.,data=trainFrame,neurons=5, epochs=1000, cores=3)  ## long time on server

preds_brnn<-predict(model_brnn,trainFrame)

########## ANN
#model_ANN<-neuralnet(Resp~.,data=trainFrame,linear.output=FALSE)
#preds_ANN<-predict(model_ANN,trainFrame)

########## Naive Bayes
## method=nnet, lm (0.0199), rpart, earth (0.07), Naive Bayes, LDA, QDA, ADA Not working
#model_LDA <-train(Resp~.,data=trainFrame,method="ada",trControl=trainControl("cv"))
#preds_LDA<-predict(model_LDA,testFrame)
