


#############################################
######### 1. GR paper linear Ridge  #########
#############################################
preds_GR<-predict(model_GR,testFrame)
GR_result<-eval_result(preds_GR)
GR_result$method <- "GR paper linear Ridge"
GR_result$drug <- drug
cat(paste(Sys.time(),"==========","1. GR paper linear Ridge Complete\n"))
#############################################
########## 2.  Random Forest
#############################################
preds_rf<-predict(model_rf,testFrame)
rf_result<-eval_result(preds_rf)
rf_result$method <- "Random Forest"
rf_result$drug <- drug
cat(paste(Sys.time(),"==========","2. Random Forest Complete\n"))

preds_ranger<-predict(model_ranger,testFrame)

#############################################
############## 3. Principle Component Regression
#############################################
preds_pcr<-predict(model_pcr,testFrame)
pcr_result<-eval_result(preds_pcr)
pcr_result$method <- "Principle Component Regression"
pcr_result$drug <- drug

cat(paste(Sys.time(),"==========","3. Principle Component Regression Complete\n"))
#############################################
##############  4. Partial Least Square
#############################################
preds_pls<-predict(model_pls,testFrame)
pls_result<-eval_result(preds_pls)
pls_result$method <- "Partial Least Square"
pls_result$drug <- drug
cat(paste(Sys.time(),"==========","4. Partial Least Square Complete\n"))
#############################################
############ 5. Ridge GLM penalty alpha=0
#############################################
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame[,-1]))
RidgeGLM_result <-eval_result(preds_ridgeglm)
RidgeGLM_result$method <- "Ridge GLM"
RidgeGLM_result$drug <- drug
cat(paste(Sys.time(),"==========","5. Ridge GLM penalty Complete\n"))

#############################################
############ 6. Lasso GLM penalty alpha=1
#############################################

preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame[,-1]))
Lasso_result_1<-eval_result(preds_Lasso_1)
Lasso_result_1$method <- "Lasso GLM"
Lasso_result_1$drug <- drug
cat(paste(Sys.time(),"==========","6. Lasso GLM penalty Complete\n"))

#############################################
############ 7. K-nearest neighbors (KNN) algorithm
#############################################
preds_knn<-predict(model_KNN,testFrame)
knn_result<-eval_result(preds_knn)
knn_result$method <- "K-nearest neighbors (KNN) algorithm"
knn_result$drug <- drug
cat(paste(Sys.time(),"==========","7. K-nearest neighbors (KNN) algorithm Complete\n"))
#############################################
############ 8. Support vector machine regression
#############################################
preds_svm<-predict(model_svm,testFrame)
svm_result<-eval_result(preds_svm)

svm_result$method <- "Support vector machine regression"
svm_result$drug <- drug

cat(paste(Sys.time(),"==========","8. Support vector machine regression Complete\n"))

#############################################
############ 10. Elastic Net
#############################################

preds_EN<-predict(model_EN,testFrame)
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
