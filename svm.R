#source("rework.R") get everything and continue from rework
library(doMC)
options(cores = 5)
registerDoMC()
set.seed(1)


eval_test<-function(preds){
  MAE <- mean(abs(testFrame$Resp - preds))
  SSE <- sum((testFrame$Resp - preds)^2)
  SST <- sum((testFrame$Resp - mean(testFrame$Resp))^2)
  SSM <- sum((preds-mean(testFrame$Resp))^2)
  R_square <- 1 - SSE / SST
  R2 <- (cor(as.numeric(preds),testFrame$Resp))^2
  RMSE = sqrt(SSE/nrow(testFrame))
  results<-list(R2=R2,RMSE=RMSE,R_Square=R_square,MAE=MAE)
  return(results)
}

for (i in 1:length(possibleDrugs2)){
drug <-  possibleDrugs2[i]
drug_data <- getDrugData(drug)

library(LiblineaR)
library(kernlab)

train_length <- round(dim(drug_data)[1]*0.8)
trainFrame <- drug_data[1:train_length,]
testFrame <- drug_data[(train_length+1):dim(drug_data)[1],]

svm <- list()
svm0<-train(Resp~.,data=trainFrame,method = 'svmLinear2',importance=TRUE)
svm1<-train(Resp~.,data=trainFrame,method = 'svmLinear3',importance=TRUE)
svm2<-train(Resp~.,data=trainFrame,method = 'svmBoundrangeString',importance=TRUE)
svm3<-train(Resp~.,data=trainFrame,method = 'svmLinear',importance=TRUE)
svm4<-train(Resp~.,data=trainFrame,method = 'svmPoly',importance=TRUE)
svm5<-train(Resp~.,data=trainFrame,method = 'svmRadial',importance=TRUE)
svm6<-train(Resp~.,data=trainFrame,method = 'svmRadialCost',importance=TRUE)
svm7<-train(Resp~.,data=trainFrame,method = 'svmRadialSigma',importance=TRUE)
svm8<-train(Resp~.,data=trainFrame,method = 'svmSpectrumString',importance=TRUE)

pred <- list()
result <-list()
  pred[[2]]<-predict(svm0,testFrame)
  pred[[1]]<-predict(svm1,testFrame)
  pred[[3]]<-predict(svm3,testFrame)
  pred[[4]]<-predict(svm4,testFrame)
  pred[[5]]<-predict(svm5,testFrame)
  pred[[6]]<-predict(svm6,testFrame)
  pred[[7]]<-predict(svm7,testFrame)

for (i in 1:7){
  cat(i)
  result[[i]]<-eval_result(unlist(pred[[i]]))
}
Result_final <- ldply(result, data.frame)
Result_final$drug <- drug

#  R2      RMSE  R_Square Adjusted_R2       MAE
#  1 0.9937746 0.1198335 0.9904886    1.000342 0.1182937
#  2 0.9937746 0.1198335 0.9904886    1.000342 0.1182937
#  3 0.9616140 0.2749934 0.9499123    1.001801 0.1690195
#  4 0.8056815 0.6631484 0.7087216    1.010473 0.4247999
#  5 0.8094000 0.6581409 0.7131040    1.010315 0.4203313
#  6 0.7275330 0.7612143 0.6162038    1.013799 0.5120520

name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output_svm_05272021/",drug,".csv", sep="",collapse = NULL)
write.csv(Result_final,name, quote=FALSE)
}
