#==========================================================================================================
#========   2.1 Validation using "80/20" cross validation in training data only ===========================
#==========================================================================================================
library(doMC)
options(cores = 35)
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

PCC_table <- data.frame()
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
  model_rf <- train(Resp~.,data=trainFrame,
                           method = "rf",
                           ntree = 50,
                           tuneGrid = data.frame(mtry = 2),
                           nodesize = 5,
                           importance = TRUE,
                           metric = "RMSE",
                           trControl = trainControl(method = "oob", seed = c(1,1)),
                           allowParallel = FALSE)
  model_ranger<-ranger(Resp~.,data=trainFrame,num.trees = 50,splitrule="variance",seed =1)
  model_pcr<-train(Resp~.,data=trainFrame,method="pcr")
  model_pls<-train(Resp~.,data=trainFrame,method="pls")
  model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))
  model_KKNN<-train(Resp~.,data=trainFrame,method="kknn",trControl=trainControl("cv",number=10))
  model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2')
  model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag')
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
preds_rf<-predict(model_rf,testFrame)
preds_ranger<-predict(model_ranger,testFrame)

preds_pcr<-predict(model_pcr,testFrame)
preds_pls<-predict(model_pls,testFrame)
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame))
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame))
preds_knn<-predict(model_KNN,testFrame)
preds_svm<-predict(model_svm,testFrame)
preds_treebag<-predict(model_treebag,testFrame)
preds_EN<-predict(model_EN,testFrame)
preds_kknn<-predict(model_KKNN,testFrame)

preds_rfinv1 <- model_rfinv1$testPred
preds_rfinv2 <- model_rfinv2$testPred
preds_rfinv3 <- model_rfinv3$testPred

preds_from_testFrame<- list(preds_GR,
            preds_rf,
            preds_ranger$predictions,
            preds_pcr,
            preds_pls,
            preds_ridgeglm,
            preds_Lasso_1,
            preds_knn,
            preds_kknn,
            preds_svm,
            preds_treebag,
            preds_EN,
            model_rfinv1$testPred,
            model_rfinv2$testPred,
            model_rfinv3$testPred)

#========================================
#  End of  ==1.0.1. Building Models==
#========================================

  for (i in 1:length(method_list)){
    PCC_table0 <- data.frame(
      drug = drug,
      method = method_list[i],
      PCC_CV28=cor(as.numeric(testResp), as.numeric(preds_from_testFrame[[i]]),method = "pearson")
    )
    PCC_table <- rbind(PCC_table,PCC_table0)
  }

  cat(paste(Sys.time(),drug," PCC finished.\n", sep="",collapse = NULL))
  file <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/PCC_CV28_eachdrug/",drug,"_PCC_CV28.csv", sep="",collapse = NULL)
  write.csv(PCC_table,file, quote=FALSE)

}


#========================================
#                 Plotting
#========================================
library("ggpubr")
library(ggplot2)
require(gridExtra)
library(cowplot)

PCC_table <- read.csv("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/39drugs_PCC_CV28.csv", header = TRUE,na.strings="NA")


filename = paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/39drugs_PCC_CV28_boxplot07162021.pdf", sep="",collapse = NULL)
pdf(filename,width=10, height=10)
ggviolin(subset(PCC_table, !is.na(PCC_CV28)), x = "method", y = "PCC_CV28",
          color = "method", ylim = c(0, 0.75),
          ylab = "PCC (Closer to 1 is better)", xlab = "method",
         repel=TRUE,draw_quantiles = 0.5,trim=FALSE,add = "jitter")+  # outlier.shape = NA
  scale_x_discrete(guide = guide_axis(angle = 20))+
  stat_summary(fun.data = function(x) data.frame(y=0.05, label =signif(median(x),digits = 3)), geom="text") +
  theme(legend.position="none")

ggboxplot(PCC_table, x = "method", y = "PCC_CV28",
          color = "method", ylim = c(0.1, 0.8),font.label = list(size = 4),
          ylab = "Pearson Cor Coef \n(Closer to 1 is better)", xlab = "method",repel=TRUE)+
          scale_x_discrete(guide = guide_axis(angle = 20))+ rremove("legend")
dev.off()
