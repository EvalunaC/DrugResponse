#source("rework.R") get everything and continue from rework

CV_methods_training <- function(drug_data,drug){
set.seed(1)
################################
####   Split the data       ####
################################
train_length <- round(dim(drug_data)[1]*0.8)
trainFrame <- drug_data[1:train_length,]
testFrame <- drug_data[(train_length+1):dim(drug_data)[1],]
################################
####   Training   Models    ####
################################

  model_GR <- linearRidge(Resp ~ ., data = trainFrame)
  model_rf <- train(Resp~.,data=trainFrame,
                           method = "rf",
                           ntree = 50,
                           tuneGrid = data.frame(mtry = 2),
                           nodesize = 5,
                           importance = TRUE,
                           metric = "RMSE",
                           trControl = trainControl(method = "oob", seed = c(1,1)),
                           allowParallel = FALSE)
  model_pcr<-train(Resp~.,data=trainFrame,method="pcr",importance=TRUE)
  model_pls<-train(Resp~.,data=trainFrame,method="pls",importance=TRUE)

  cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
  (best_lam_0 <- cv_output_0$lambda.min)
  model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)
  cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
  (best_lam_1 <- cv_output_1$lambda.min) ### 0.1444154
  model_Lasso_1<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1, lambda=best_lam_1)
  model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
  model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2',importance=TRUE)
  model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))


methods_model<- list(model_GR,model_rf,model_pcr,model_pls,model_KNN,model_svm,model_EN,model_ridgeglm,model_Lasso_1)
methods <- c("GR paper linear Ridge",
            "Random Forest",
            "Principle Component Regression",
            "Partial Least Square",
            "KNN",
            "SVM",
            "Elastic Net",
            "Ridge GLM",
            "Lasso GLM")



################################
####  Testing Predictions   ####
################################
cv28_pcc <- c()
#pctest_p <- c()
for (i in 1:7){
  cat(i)
  preds<-predict(methods_model[i],testFrame)
  cv28_pcc[i] <- cor(as.numeric(testFrame$Resp), as.numeric(unlist(preds)),method = "pearson")
}
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame[,-1]))
cv28_pcc[8] <- cor(as.numeric(testFrame$Resp), as.numeric(unlist(preds_ridgeglm)),method = "pearson")
preds_Lasso <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame[,-1]))
cv28_pcc[9] <- cor(as.numeric(testFrame$Resp), as.numeric(unlist(preds_Lasso)),method = "pearson")

result <- t(rbind(drug,methods,cv28_pcc))
name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output_cv28pcc_05132021/",drug,".csv", sep="",collapse = NULL)
write.csv(result,name, quote=FALSE)
}

#run the function
for (i in 2:length(possibleDrugs2)){
  drug_data <- getDrugData(possibleDrugs2[i])
  CV_methods_training(drug_data,possibleDrugs2[i])
  }

  #====================================================================================================
  # showing result
  #====================================================================================================
  result_list <- list.files(path="/extraspace/ychen42/Drug_Response/yiqings_work/Output_cv28pcc_05132021/")
  length(result_list)
  result0 <- data.frame()

  length(result_list)
  for (i in 1:length(result_list)){
    file <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output_cv28pcc_05132021/",result_list[i], sep="",collapse = NULL)
    result_file <- read.csv(file, header = TRUE,na.strings="NA")
    result0 <- rbind(result0,result_file)
  }
