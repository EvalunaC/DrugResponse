

library(ggplot2)

library(tidyr)
library(pROC)
library(dplyr )

neural.network <- read.csv("/home/ychen42/Lapatinib_y_hat.csv",header=FALSE)
colnames(neural.network) <- "neural.network"
y <- read.csv("/extraspace/ychen42/Drug_Response/Data/DLprediction/Lapatinib_y.csv",header=FALSE)
y_hat <- read.csv("/extraspace/ychen42/Drug_Response/Data/DLprediction/Lapatinib_y_hat.csv",header=FALSE)
pred <- read.csv("/extraspace/ychen42/Drug_Response/Data/DLprediction/Lapatinib_y_pred.csv",header=FALSE)
data <- data.frame(y,y_hat)

data$y_cate <- ifelse(data$V1>=mean(data$V1),1,0)
data$yhat_cate <- ifelse(data$V1.1>=mean(data$V1.1),1,0)
(roc.pred2 <- roc(data$y_cate, data$yhat_cate, percent = TRUE, main = "Smoothing"))


pdf("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Lapatinib_scatter.pdf")
plot(data$V1, data$V1.1, pch = 19, col = "lightblue",main = "Lapatinib prediction vs. true response")
# Regression line
abline(lm(data$V1.1 ~ data$V1), col = "red", lwd = 3)
# Pearson correlation
text(paste("Correlation:", round(cor(data$V1, data$V1.1), 2)), x = 25, y = 95)
plot.roc(data$y_cate, data$yhat_cate, percent = TRUE, main = "Categorized Lapatinib prediction vs. true response ROC curves", print.auc = TRUE)
dev.off()

abline(lm(data$V1.1 ~ data$V1), col = "red", lwd = 3)


preds_from_trainFrame<- data.frame(trainFrame$Resp,
            preds_GR,
            preds_rf,
            preds_ranger$predictions,
            preds_pcr,
            preds_pls,
            preds_ridgeglm,
            preds_Lasso_1,
            preds_knn,
            preds_kknn,
            preds_svm,
            preds_EN,
            model_rfinv1$testPred,
            model_rfinv2$testPred,
            model_rfinv3$testPred,
          neural.network)
write.csv(preds_from_trainFrame,"/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Lapatinib_prediction.csv",quote = FALSE)
predict <- tidyr::gather(preds_from_trainFrame, method, predict, preds_GR:model_rfinv3$testPred)






#########plot#############
predict <- read.csv("C:/Users/EYC/Downloads/Lapatinib_prediction.csv",header=T)

predict <- predict %>%
  rename(GLMRidge=X1 , GLMLasso=X1.1)

predict_melt <- tidyr::gather(predict, method, predict, preds_GR:neural.network)


ggplot(predict_melt, aes(x = trainFrame.Resp, y = predict, colour = method)) + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                                                                           linetype="dashed", size=1)+
  geom_point() + xlab('response level')+facet_wrap( ~ method, nrow = 3)

########## AUC ROC ###############
predict_melt$predict_cate <- ifelse(predict_melt$predict>=mean(predict_melt$predict),1,0)
predict_melt$Resp_cate <- ifelse(predict_melt$trainFrame.Resp>=mean(predict_melt$trainFrame.Resp),1,0)
(roc.pred2 <- roc(data$y_cate, data$yhat_cate, percent = TRUE, main = "Smoothing"))


pdf("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Lapatinib_scatter.pdf")
plot(data$V1, data$V1.1, pch = 19, col = "lightblue",main = "Lapatinib prediction vs. true response")
# Regression line
abline(lm(data$V1.1 ~ data$V1), col = "red", lwd = 3)
# Pearson correlation
text(paste("Correlation:", round(cor(data$V1, data$V1.1), 2)), x = 25, y = 95)
plot.roc(data$y_cate, data$yhat_cate, percent = TRUE, main = "Categorized Lapatinib prediction vs. true response ROC curves", print.auc = TRUE)
dev.off()



library(pROC)
methods <- unique(predict_melt$method)
roc <- list()
data <- predict_melt[predict_melt$method==methods[1],]
percent  <- c()
plot.roc(data$Resp_cate, data$predict_cate, percent = F, main = "ROC curves", add =  FALSE, asp = NA, print.auc = TRUE)
for (i in 2:length(methods)){
  data <- predict_melt[predict_melt$method==methods[i],]
  roc.model <- roc(data$Resp_cate, data$predict_cate, data = data, percent = F, main = "Smoothing")
  roc <- rbind(roc,roc.model)
  lines(roc.model, type = "l", lty = i, col = i)
  percent[i]<- roc.model$percent
  
}

legend(list(x = 0.4,y = 0.4), 
       legend = method_list, 
       col = 1:15,
       lty = 1:15,)



lines( type = "l", lty = 2, col = "grey35")
method_list <- c("GR paper linear Ridge",
                 "Random Forest",
                 "Random Forest (Ranger)",
                 "Principle Component Regression",
                 "Partial Least Square",
                 "Ridge GLM",
                 "Lasso GLM",
                 "K-nearest neighbors ",
                 "Weighted K-nearest neighbors ",
                 "Support Vector Machine",
                 #              "Treebag (bootstrap aggregating)",
                 "Elastic Net",
                 "RF+Lasso2019(out-of-bag)",
                 "RF+Lasso2019(split-conformal)",
                 "RF+Lasso2019(quantile regression forest)",
                 "Neural Network"
)


data <- predict_melt[predict_melt$method=="neural.network",]
roc.pred1 <- roc(data$Resp_cate, data$predict_cate, data = data, percent = TRUE, main = "Smoothing")

roc.pred3 <- roc(predict$IAC, predict$phat3, percent = TRUE, main = "Smoothing")
roc.pred4 <- roc(predict$IAC, predict$phat4, percent = TRUE, main = "Smoothing")
roc.pred5 <- roc(predict$IAC, predict$phat5, percent = TRUE, main = "Smoothing")

plot.roc(predict$IAC, predict$phat1, percent = TRUE, main = "ROC curves", add =  FALSE, asp = NA, print.auc = TRUE)

lines(roc.pred2, type = "l", lty = 2, col = "grey35")
lines(roc.pred3, type = "l", lty = 3, col = "grey48")
lines(roc.pred4, type = "l", lty = 4, col = "grey61")
lines(roc.pred5, type = "l", lty = 1, pch = 24, col = "grey76")

legend("bottomright", 
       legend = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"), 
       col = c("black", "grey35", "grey48", "grey61", "grey76"),
       lty = c(1, 2, 3, 4, 1))







predict <- data.frame()
for (i in 1:length(method_list)){
    predict0 <- data.frame(method = method_list[i], predict = preds_from_trainFrame[[i]])
    predict <- rbind(predict,predict0)
}


gapminder_co2 %>%
  ggplot(aes(x=gdpPercap,y=co2)) +
  geom_point(alpha=0.5, aes(color=continent,size=pop)) +
  labs(x="GDP per capita", y= "CO2 Emission per person (in tonnes)",
       title="CO2 emission per person vs GDP per capita") +
  scale_y_log10()+
  scale_x_log10()+
  geom_smooth(method=lm,se=FALSE) +
  ggrepel::geom_text_repel(data = gapminder_co2 %>%
                             filter(gdpPercap>12000  | gdpPercap < 1000) %>%
                             sample_n(20),
                           aes(label = country))

Stat_table <- data.frame()
for (i in 1:length(method_list)){
    Stat_table0 <- eval_result_train(preds=preds_from_trainFrame[[i]])
    Stat_table0$method <- method_list[i]
    Stat_table <- rbind(Stat_table,Stat_table0)
}





setwd("/extraspace/ychen42/Drug_Response/Data/")

GDSC2<-read.csv("GDSC2_IC50_matrix.csv", header=T,row.names=1)

#GDSC2<-read.csv("GDSC2_matrix.csv", header=T,row.names=1)
possibleDrugs2<-rownames(GDSC2) ##192
drug_list <-possibleDrugs2[1:20]


#source("Creating_drug_data.R")
memory.limit(9999999999)


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

library(ranger)
library(kknn)
library(rfinterval)

#drug <- "Lapatinib"
#drug_data <- getDrugData("Lapatinib")
#drug_data <- getDrugData(possibleDrugs2[i])
#trainFrame <- drug_data[[1]]
#testFrame <- drug_data[[2]]

# To run: build_stateval(drug=drug,trainFrame=drug_data[[1]],testFrame=drug_data[[2]])

method_list <- c("GR paper linear Ridge",
              "Random Forest",
              "Random Forest (Ranger)",
              "Principle Component Regression",
              "Partial Least Square",
              "Ridge GLM",
              "Lasso GLM",
              "K-nearest neighbors ",
              "Weighted K-nearest neighbors ",
              "Support Vector Machine",
#              "Treebag (bootstrap aggregating)",
              "Elastic Net",
              "RF+Lasso2019(out-of-bag)",
              "RF+Lasso2019(split-conformal)",
              "RF+Lasso2019(quantile regression forest)"
              )

build_stateval <- function(drug,trainFrame,testFrame=trainFrame[,-1],outFolder){ # trainFrame[,-1] is default test set

  #==================================================
  #========   1.0.0. Validation Function
  #==================================================
  eval_result_train<-function(preds){
    MAE <- mean(abs(trainFrame$Resp - preds))
    SSE <- sum((trainFrame$Resp - preds)^2)
    SST <- sum((trainFrame$Resp - mean(trainFrame$Resp))^2)
    SSM <- sum((preds-mean(trainFrame$Resp))^2)
    R2_simp <- 1 - SSE / SST
    R2_corr <- (cor(as.numeric(preds),trainFrame$Resp))^2
    RMSE = sqrt(SSE/nrow(trainFrame))
    results<-data.frame(RMSE=RMSE,MAE=MAE,R2_simp=R2_simp,R2_corr=R2_corr)
    return(results)
  }

  eval_result_test<-function(preds){
      MAE <- mean(abs(trainFrame$Resp - preds))
      SSE <- sum((trainFrame$Resp - preds)^2)
      SST <- sum((trainFrame$Resp - mean(trainFrame$Resp))^2)
      SSM <- sum((preds-mean(trainFrame$Resp))^2)
      R2_simp <- 1 - SSE / SST
      R2_corr <- (cor(as.numeric(preds),trainFrame$Resp))^2
      RMSE = sqrt(SSE/nrow(trainFrame))
      results<-data.frame(RMSE=RMSE,MAE=MAE,R2_simp=R2_simp,R2_corr=R2_corr)
      return(results)
    }
  #==================================================
  #========    1.0.1. Building Models
  #==================================================
set.seed(1)

cat(paste(Sys.time(),"==========",drug,":\n"))

cat(paste(Sys.time(),"==========","GR paper linear Ridge Start...\n"))
model_GR <- linearRidge(Resp ~ ., data = trainFrame)

cat(paste(Sys.time(),"==========","Random Forest Start...\n"))
model_rf <- train(Resp~.,data=trainFrame,
                         method = "rf",
                         ntree = 50,
                         tuneGrid = data.frame(mtry = 2),
                         nodesize = 5,
                         importance = TRUE,
                         metric = "RMSE",
                         trControl = trainControl(method = "oob", seed = c(1,1)),
                         allowParallel = FALSE)

cat(paste(Sys.time(),"==========","Random Forest Ranger Start...\n"))
#model_ranger<-train(Resp~.,data=trainFrame,num.trees = 50,method="ranger",trControl = trainControl(method = "oob", seed = list(c(1,1,1,1,1,1),c(1,1,1,1,1,1))))
model_ranger<-ranger(Resp~.,data=trainFrame,num.trees = 50,splitrule="variance",seed =1)

cat(paste(Sys.time(),"==========","Principle Component Regression Start...\n"))
model_pcr<-train(Resp~.,data=trainFrame,method="pcr")#,importance=TRUE)

cat(paste(Sys.time(),"==========","Partial Least Square Start...\n"))
model_pls<-train(Resp~.,data=trainFrame,method="pls")#,importance=TRUE)

cat(paste(Sys.time(),"==========","K-nearest neighbors (KNN) algorithm Start...\n"))
model_KNN<-train(Resp~.,data=trainFrame,method="knn",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","Weighted K-nearest neighbors (KNN) algorithm Start...\n"))
model_KKNN<-train(Resp~.,data=trainFrame,method="kknn",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","Support vector machine regression Start...\n"))
model_svm<-train(Resp~.,data=trainFrame,method = 'svmLinear2')#,importance=TRUE)

#cat(paste(Sys.time(),"==========","Treebag (bootstrap aggregating) algorithm Start...\n"))
#model_treebag<-train(Resp~.,data=trainFrame,method = 'treebag')

cat(paste(Sys.time(),"==========","Elastic Net Regression Start...\n"))
model_EN<-train(Resp~.,data=trainFrame,method="glmnet",trControl=trainControl("cv",number=10))

cat(paste(Sys.time(),"==========","Ridge GLM penalty Start...\n"))
cv_output_0 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0,type.measure="mse",nfolds=10)
(best_lam_0 <- cv_output_0$lambda.min) ###104.2796    Evry time different. 49.54127, 47.28955, 114.4468
model_ridgeglm<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=0, lambda=best_lam_0)

cat(paste(Sys.time(),"==========","Lasso GLM penalty Start...\n"))
cv_output_1 <- cv.glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1,type.measure="mse",nfolds=10)
(best_lam_1 <- cv_output_1$lambda.min) ### 0.1444154
model_Lasso_1<- glmnet(as.matrix(trainFrame[,-1]),as.matrix(trainFrame$Resp),alpha=1, lambda=best_lam_1)

#### Regression-Enhanced Random Forests, "Predictive Inference for Random Forests" RF+Lasso
cat(paste(Sys.time(),"==========","Predictive Inference for Random Forests (oob) Start...\n"))

model_rfinv1<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
                         method = "oob", alpha = 0.05,
                         symmetry = TRUE)

cat(paste(Sys.time(),"==========","Predictive Inference for Random Forests (split-conformal) Start...\n"))
model_rfinv2<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
                         method = "split-conformal", alpha = 0.05,
                         seed = 1)
cat(paste(Sys.time(),"==========","Predictive Inference for Random Forests (quantreg) Start...\n"))
model_rfinv3<- rfinterval(Resp~.,train_data=trainFrame, test_data = testFrame,
                         method = "quantreg", alpha = 0.05)

  #==================================================
  #===       1.0.2. Prediction & Evaluation
  #==================================================


#===========  1.0.2.1 Using trainFrame   =========================================

preds_GR<-predict(model_GR,trainFrame[,-1])
preds_rf<-predict(model_rf,trainFrame[,-1])
preds_ranger<-predict(model_ranger,trainFrame[,-1])

preds_pcr<-predict(model_pcr,trainFrame[,-1])
preds_pls<-predict(model_pls,trainFrame[,-1])
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(trainFrame[,-1]))
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(trainFrame[,-1]))
preds_knn<-predict(model_KNN,trainFrame[,-1])
preds_svm<-predict(model_svm,trainFrame[,-1])
#preds_treebag<-predict(model_treebag,trainFrame[,-1])
preds_EN<-predict(model_EN,trainFrame[,-1])
preds_kknn<-predict(model_KKNN,trainFrame[,-1])

preds_rfinv1 <- model_rfinv1$testPred
preds_rfinv2 <- model_rfinv2$testPred
preds_rfinv3 <- model_rfinv3$testPred

preds_from_trainFrame<- list(preds_GR,
            preds_rf,
            preds_ranger$predictions,
            preds_pcr,
            preds_pls,
            preds_ridgeglm,
            preds_Lasso_1,
            preds_knn,
            preds_kknn,
            preds_svm,
#            preds_treebag,
            preds_EN,
            model_rfinv1$testPred,
            model_rfinv2$testPred,
            model_rfinv3$testPred)

Stat_table <- data.frame()
for (i in 1:length(method_list)){
    Stat_table0 <- eval_result_train(preds=preds_from_trainFrame[[i]])
    Stat_table0$method <- method_list[i]
    Stat_table <- rbind(Stat_table,Stat_table0)
}
Stat_table$drug <- drug
df <- Stat_table[,c("method","drug","RMSE","MAE","R2_simp","R2_corr")]
name <- paste(outFolder,drug,"_4Stats.csv", sep="",collapse = NULL)
write.csv(df,name, quote=FALSE)

cat(paste(Sys.time(),"Stat_table saved as ",outFolder,"/..._4Stats.csv \n", sep="",collapse = NULL))
}


###===========  1.0.2.2 Potential: Using testFrame   =========================================
preds_GR<-predict(model_GR,testFrame)
preds_rf<-predict(model_rf,testFrame)
preds_ranger<-predict(model_ranger,testFrame)

preds_pcr<-predict(model_pcr,testFrame)
preds_pls<-predict(model_pls,testFrame)
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame))
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame))
preds_knn<-predict(model_KNN,testFrame)
preds_svm<-predict(model_svm,testFrame)
#preds_treebag<-predict(model_treebag,testFrame)
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

  Stat_table <- data.frame()
  for (i in 1:length(method_list)){
    Stat_table0 <- eval_result_train(preds=preds_from_testFrame[[i]])
    Stat_table0$method <- method_list[i]
    Stat_table <- rbind(Stat_table,Stat_table0)
  }
  Stat_table$drug <- drug
  #  Stat_table$test_data <- deparse(substitute(x))
  #df <- Stat_table[,c("method","drug",,"test_data","RMSE","MAE","R2_simp","R2_corr")]
  df <- Stat_table[,c("method","drug","RMSE","MAE","R2_simp","R2_corr")]
  name <- paste(outFolder,drug,"_4Stats.csv", sep="",collapse = NULL)
  write.csv(df,name, quote=FALSE)

  cat(paste(Sys.time(),"Stat_table saved as ",outFolder,"/..._4Stats.csv \n", sep="",collapse = NULL))

#filename = paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Models/",drug,"_preds_from_testFrame.RData", sep="",collapse = NULL)
#save(preds_from_testFrame, method_list, file=filename)

#cat(paste(Sys.time(),"testFrame pred saved as /extraspace/ychen42/Drug_Response/Own_update2.0/Models/..._preds_from_testFrame.RData\n"))

#return(list(preds_from_trainFrame,preds_from_testFrame,Stat_table))
}

#====================  End of function "build_stateval" ===============

#==================================================
#===   Potential:    1.0.5. Saving models (use in function)
#==================================================
#setwd("/extraspace/ychen42/Drug_Response/Own_update2.0/Models/")
#
#saveRDS(model_GR,paste(drug,"_model_GR.rds", sep="",collapse = NULL))
#saveRDS(model_rf,paste(drug,"_model_rf.rds", sep="",collapse = NULL))
#saveRDS(model_ranger,paste(drug,"_model_ranger.rds", sep="",collapse = NULL))
#saveRDS(model_pcr,paste(drug,"_model_pcr.rds", sep="",collapse = NULL))
#saveRDS(model_pls,paste(drug,"_model_pls.rds", sep="",collapse = NULL))
#saveRDS(model_ridgeglm,paste(drug,"_model_ridgeglm.rds", sep="",collapse = NULL))
#saveRDS(model_Lasso_1,paste(drug,"_model_Lasso_1.rds", sep="",collapse = NULL))
#saveRDS(model_KNN,paste(drug,"_model_KNN.rds", sep="",collapse = NULL))
#saveRDS(model_svm,paste(drug,"_model_svm.rds", sep="",collapse = NULL))
#saveRDS(model_treebag,paste(drug,"_model_treebag.rds", sep="",collapse = NULL))
#saveRDS(model_EN,paste(drug,"_model_EN.rds", sep="",collapse = NULL))
#saveRDS(model_KKNN,paste(drug,"_model_KKNN.rds", sep="",collapse = NULL))
#saveRDS(model_rfinv1,paste(drug,"_model_rfinv1.rds", sep="",collapse = NULL))
#saveRDS(model_rfinv2,paste(drug,"_model_rfinv2.rds", sep="",collapse = NULL))
#saveRDS(model_rfinv3,paste(drug,"_model_rfinv3.rds", sep="",collapse = NULL))

#To read: super_model <- readRDS("./final_model.rds")
