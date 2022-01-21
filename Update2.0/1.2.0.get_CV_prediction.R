
method_list <- c("GR paper linear Ridge",
#              "Random Forest",
              "Random Forest (Ranger)",
              "Principle Component Regression",
              "Partial Least Square",
              "Ridge GLM",
              "Lasso GLM",
              "K-nearest neighbors ",
#              "Weighted K-nearest neighbors ",
              "Support Vector Machine",
#              "Treebag (bootstrap aggregating)",
              "Elastic Net",
              "RF+Lasso2019(out-of-bag)",
              "RF+Lasso2019(split-conformal)",
              "RF+Lasso2019(quantile regression forest)"
              )
library(doMC)
options(cores = 45)
registerDoMC()

#source("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/0.1.Creating_drug_data.R")
#source("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/1.0.build_models_stateva.R")

set.seed(1)

#for (i in 1:length(possibleDrugs2)){
#  drug_data <- getDrugData(possibleDrugs2[i])
#  trainFrame <-
#  drug<-possibleDrugs2[i]
#drug <- "Lapatinib"
#drug_data <- getDrugData("Lapatinib")[[1]] #trainFrame only
#load("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/Lapatinib_trainFrame.RData")
#load("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/Lapatinib_testFrame.RData")

prediction <- data.frame()
for (i in 2:length(possibleDrugs2)){
  drug <- possibleDrugs2[i]
  cat(paste(Sys.time()," Reading ",drug,"...\n", sep="",collapse = NULL))
file <- paste("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/",possibleDrugs2[i],"_trainFrame.RData", sep="",collapse = NULL)
load(file)
drug_data <- trainFrame
  #========================================
  #=====    Split the data     ============
  #========================================
  train_length <- round(dim(drug_data)[1]*0.8)
  trainFrame <- drug_data[1:train_length,]
  testFrame <- drug_data[(train_length+1):dim(drug_data)[1],-1]
  testResp <- drug_data[(train_length+1):dim(drug_data)[1],1]
  #========================================
  #  Then use  ==1.0.1. Building Models==
  #========================================
  set.seed(1)

  model_GR <- linearRidge(Resp ~ ., data = trainFrame)
  model_ranger<-ranger(Resp~.,data=trainFrame,num.trees = 50,splitrule="variance",seed =1)
  model_pcr<-train(Resp~.,data=trainFrame,method="pcr")
  model_pls<-train(Resp~.,data=trainFrame,method="pls")
  model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
  model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2')
  model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))
  cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
  (best_lam_0 <- cv_output_0$lambda.min)
  model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)
  cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
  (best_lam_1 <- cv_output_1$lambda.min)
  model_Lasso_1<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1, lambda=best_lam_1)
  model_rfinv1<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
                           method = "oob", alpha = 0.05,
                           symmetry = TRUE)
  model_rfinv2<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
                           method = "split-conformal", alpha = 0.05,
                           seed = 1)
  model_rfinv3<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
                           method = "quantreg", alpha = 0.05)

preds_GR<-predict(model_GR,testFrame)
preds_ranger<-predict(model_ranger,testFrame)
preds_pcr<-predict(model_pcr,testFrame)
preds_pls<-predict(model_pls,testFrame)
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame))
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame))
preds_knn<-predict(model_KNN,testFrame)
preds_svm<-predict(model_svm,testFrame)
preds_EN<-predict(model_EN,testFrame)
preds_rfinv1 <- model_rfinv1$testPred
preds_rfinv2 <- model_rfinv2$testPred
preds_rfinv3 <- model_rfinv3$testPred

preds_from_testFrame<- data.frame(preds_GR,
            preds_ranger$predictions,
            preds_pcr,
            preds_pls,
            preds_ridgeglm,
            preds_Lasso_1,
            preds_knn,
            preds_svm,
            preds_EN,
            model_rfinv1$testPred,
            model_rfinv2$testPred,
            model_rfinv3$testPred)

  colnames(preds_from_testFrame) <- c("GR paper linear Ridge",
                "Random Forest (Ranger)",
                "Principle Component Regression",
                "Partial Least Square",
                "Ridge GLM",
                "Lasso GLM",
                "K-nearest neighbors ",
  #              "Weighted K-nearest neighbors ",
                "Support Vector Machine",
  #              "Treebag (bootstrap aggregating)",
                "Elastic Net",
                "RF+Lasso2019(out-of-bag)",
                "RF+Lasso2019(split-conformal)",
                "RF+Lasso2019(quantile regression forest)"
                )
library(reshape2)
preds_from_testFrame$id <- rownames(preds_from_testFrame)
preds_long <- melt(preds_from_testFrame,id="id")
preds_long$drug_name <- drug

prediction <- rbind(prediction,preds_long)
cat(paste(Sys.time(),"==========",drug,"\n"))
file <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Output/192CV_prediction/",drug,".csv", sep="",collapse = NULL)
write.csv(prediction,file, quote=FALSE)
}
