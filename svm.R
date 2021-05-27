#source("rework.R") get everything and continue from rework
drug_data <- getDrugData("Lapatinib")
trainFrame <- drug_data


library(LiblineaR)
library(kernlab)


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
for (i in 3:7){
  cat(i)
  pred[[i]]<-predict(unlist(svm[[i]]),trainFrame)
  result[i]<-eval_result(unlist(pred[i]))
}
result


preds<-predict(svm1,trainFrame)
result[[1]]<-eval_result(preds)
