set.seed(1000)
memory.limit(9999999999)

load("C:/Users/qiangli/Desktop/Drug Response Project/DR_Pred/trainFrame.RData")

library(caret)
library(pRRophetic)
library(ridge)
library(MASS)
library(dplyr)
library(caret)
library(glmnet)

### RMSE and R2 function
eval_result<-function(preds){
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


######### GR paper linear Ridge
model_GR <- linearRidge(Resp ~ ., data = trainFrame)
preds_GR<-predict(model_GR,trainFrame)

GR_result<-eval_result(preds_GR) ## RMSE=0.639 R2=0.729 R2_adj=1.01 MAE=0.499 AIC=

result_model_GR_nPC<-model_GR$df[model_GR$chosen.nPCs,]  #97
#    model  variance  residual
#195.24401  97.75603 292.73200

summary_GRsummary(model_GR)$summaries$summary97


 ####### Ridge model by lm.ridge ????? coef is extremely samll

model_Ridge<-lm.ridge(Resp ~ ., data = trainFrame)
# Ridge_preds<-predict(model_Ridge,trainFrame)  ## Error
preds_Ridge<-as.matrix(cbind(const=1,trainFrame))[,-2] %*% coef(model_Ridge)

Ridge_result<-eval_result(preds_Ridge)


############## Principle Component Regression
library(pls)

model_pcr<-pcr(Resp~.,data=trainFrame,ncomp=3, validation = "CV", jackknife = TRUE)
#jack.test(model_pcr, ncomp = 3)
preds_pcr<-predict(model_pcr,trainFrame,type="response")

# eval_result(preds_pcr)  #generate wrong number R2<0
rmsep_pcr <- sqrt(mean((trainFrame$Resp - preds_pcr)^2)) ##1.1897

#cor(as.numeric(preds_pcr),brcaLap,method="spearman")
obsfit <- predplot(model_pcr, which = "validation")
Residuals <- obsfit[,1] - obsfit[,2]


##############   Partial Least Square

model_plsr<-plsr(Resp~.,data=trainFrame, ncomp=3, validation="CV",jackknife=TRUE)
#jack.test(model_plsr,ncomp=3)
preds_plsr<-predict(model_plsr,trainFrame,type="response")

# eval_result(preds_plsr)  #generate wrong number R2<0
rmsep_plsr <- sqrt(mean((trainFrame$Resp - preds_plsr)^2)) ##0.9489


############ Lasso
 #without penalty: alpha=0 (Risge Penalty)
cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
(best_lam_0 <- cv_output_0$lambda.min) ###104.2796
model_Lasso_0<- glmnet(trainFrame[,-1],trainFrame$Resp,alpha=0, lambda=best_lam_0)
preds_Lasso_0 <- predict(model_Lasso_0, s = best_lam_0, newx=as.matrix(trainFrame[,-1]))

Lasso_result_0<-eval_result(preds_Lasso_0)
AIC_BIC_Lasso_0<-AIC_BIC_glmnet(model_Lasso_0)

#without penalty: alpha=1
cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
(best_lam_1 <- cv_output_1$lambda.min) ### 0.1444154
model_Lasso_1<- glmnet(trainFrame[,-1],trainFrame$Resp,alpha=1, lambda=best_lam_1)
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(trainFrame[,-1]))

Lasso_result_1<-eval_result(preds_Lasso_1)
AIC_BIC_Lasso_1<-AIC_BIC_glmnet(model_Lasso_1)


########## Random Forest
library(randomForest)

model_RF<- randomForest(Resp~.,data=trainFrame, importance = TRUE)
preds_RF<-predict(model_RF,trainFrame)


RF_result<-eval_result(preds_RF)
cr_RF<-rfcv(trainFrame[,-1],trainFrame$Resp,cv.fold=10)


###### SVM
library(e1071)

model_SVM<-svm(Resp~.,data=trainFrame)
preds_SVM<-predict(model_SVM,trainFrame)

SVM_result<-eval_result(preds_SVM)

############# Treebag
model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag',trControl=trainControl("cv",number=10))
preds_treebag<-predict(model_treebag,trainFrame)

result_model_treebag<-model_treebag$results ## RMSE=1.09846, R2=0.212 MAE=0.8788
treebag_result<-eval_result(preds_treebag)

#AIC_BIC_treebag<-AIC_BIC_train(model_treebag)




############# Elastic Net
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))
(model_EN$bestTune)
preds_EN<-predict(model_EN,trainFrame)

EN_result<-eval_result(preds_EN) ## RMSE=0.98 R2=0.363 R2_adj=1.02 MAE=0.7977 AIC=

model_EN$results$RMSE[6] ## 1.1026
model_EN$results$Rsquared[6] ##0.2227984

AIC_BIC_EN<-AIC_BIC_train(model_EN)

######## KNN
model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
preds_KNN<-predict(model_KNN,trainFrame)

KNN_result<-eval_result(preds_KNN)
result_model_KNN<-model_KNN$results

#AIC_BIC_KNN<-AIC_BIC_train(model_KNN)

#################### the following methods generate errors or taking forever

########### Robust Statistics  ## can not run
#robfit <- glmrob(Resp~.,data=trainFrame,family=gaussian,method="Mqle")
#plot(robfit)
#mean(robfit$residuals^2)
#summary(robfit)

############## Quantile Regression Neural Network ## takes too long time
# library(qrnn)
#model_qrnn<-train(Resp~.,data=trainFrame,method = 'qrnn',trControl=trainControl("cv",number=10))
#preds_qrnn<-predict(model_qrnn,trainFrame)

#model_qrnn <- qrnn.fit(x=as.matrix(trainFrame[,2:13543]),y=as.matrix(trainFrame$Resp),data=trainFrame,method="adam",n.hidden=2, n.trials=1, tau=0.2678857,n.errors.max=1000,iter.max=50,bag=TRUE) ## Error

############## ???? Robust Linear Model

#model_rlm<-train(Resp~.,data=trainFrame,method = 'rlm',trControl=trainControl("cv",number=10)) ## Error
#preds_rlm<-predict(model_rlm,testFrame)
#preds_rlm<- preds_rlm^(1/0.8944584)
#preds_rlm <- preds_rlm - 2.944905

# model_rlm<-rlm(Resp~.,data=trainFrame,psi = psi.bisquare) ##Error: x is singular


#############  Quantile Regression - conquer method

#library(quantreg)
# model_rq <- rq(Resp~.,data=trainFrame, method = "conquer") the number of conlumns of x can not exceed the number of rows x

#library(conquer)
#model_conquer <- conquer(X=as.matrix(trainFrame[,-1]),Y=trainFrame$Resp,checkSing = TRUE) #same with above


############# Bayesian Regularized Neural Networks
#library(brnn)
#library(parallel)
#model_brnn<-brnn(Resp~.,data=trainFrame,neurons=5, epochs=1000, cores=3)  ## long time on server

#preds_brnn<-predict(model_brnn,trainFrame)

########## ANN
#model_ANN<-neuralnet(Resp~.,data=trainFrame,linear.output=FALSE)
#preds_ANN<-predict(model_ANN,trainFrame)

########## Naive Bayes
## method=nnet, lm (0.0199), rpart, earth (0.07), Naive Bayes, LDA, QDA, ADA Not working
#model_LDA <-train(Resp~.,data=trainFrame,method="ada",trControl=trainControl("cv"))
#preds_LDA<-predict(model_LDA,testFrame)
