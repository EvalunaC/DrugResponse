setwd("/extraspace/ychen42/Drug_Response/Data/")
library(doMC)
options(cores = 40)
registerDoMC()


library(TCGAbiolinks)
library(SummarizedExperiment)

library(survival)
library(survminer)
library(gridExtra)
library(stringr)

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

possibleDrugs2 <- read.csv("/extraspace/ychen42/Drug_Response/Data/192Drug_list.csv")


project <- c("TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
"TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM")

surv_result <- data.frame()

#TCGA-LGG:
for (l in 3:length(project)){
#cate <- "TCGA-LGG"
cate <- project[l]
clin <- GDCquery_clinic(cate, "clinical")

query <- GDCquery(project = cate,
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab",
                legacy=TRUE)

GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
names(clinical.drug)
drug <- data.frame(submitter_id= clinical.drug$bcr_patient_barcode,drug_name = clinical.drug$drug_name)

clinic_data <- merge(clin,drug,by="submitter_id")
dim(clinic_data)
a <- table(clinic_data$drug_name)
drug_names <- intersect(names(a[a>200]),as.vector(possibleDrugs2$x))#More than 50 pts.
print(a[a>200])
print(drug_names)
#TCGA-ACC     [1] 18 73     pass
#TCGA-BLCA    [1] 303  76

for (i in 2:length(drug_names)){

  drug_data <- clinic_data[clinic_data$drug_name==drug_names[i],]
  trainFile <- paste("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/",drug_names[i],"_trainFrame.RData", sep="",collapse = NULL)
  testFile <- paste("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/",drug_names[i],"_testFrame.RData", sep="",collapse = NULL)

  load(trainFile)
  load(testFile)

  #===========================================================================
  #===========================================================================
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

#===========================================================================
#===========================================================================
preds_bi <- list()
  for (j in 1:15){
    preds_bi[[j]] <- ifelse(preds_from_testFrame[[j]] > mean(preds_from_testFrame[[j]]), "HR", "LR")
  }

  names(preds_from_testFrame) <- seq(1,15)
#  sampleNames <- colnames(tpmDatMat_bc_tpm_logged)
#  theTumorSamples <- which(substring(sampleNames, 14, 16) == "01A") # identify the tumor samples, tumor samples annotated as "01" by TCGA, normal samples as "10".
  newNames <- gsub(".", "-", substring(colnames(tpmDatMat_bc_tpm_logged), 1, 12), fixed=T)

  pred_data <- data.frame(submitter_id=newNames,
                          sampleNames=colnames(tpmDatMat_bc_tpm_logged),
                          tumor_bi=ifelse(substring(sampleNames, 14, 16) == "01A",T,F),# identify the tumor samples, tumor samples annotated as "01" by TCGA, normal samples as "10".
                          preds=preds_from_testFrame
  )
  pred_data$preds.6 <- pred_data$preds.1.1
  pred_data$preds.7 <- pred_data$preds.1.2

  pred_tumor <- pred_data[pred_data$tumor_bi==T,]# 9744    5 # Only include the tumor samples in this analysis. Results on normal samples are meaningless.
  pred_merged <- merge(pred_tumor,drug_data,by="submitter_id") #[1] 194  73
pred_merged$preds.bi.1 <- ifelse(pred_merged$preds.1 > mean(pred_merged$preds.1), "HR", "LR")
pred_merged$preds.bi.2 <- ifelse(pred_merged$preds.2 > mean(pred_merged$preds.2), "HR", "LR")
pred_merged$preds.bi.3 <- ifelse(pred_merged$preds.3 > mean(pred_merged$preds.3), "HR", "LR")
pred_merged$preds.bi.4 <- ifelse(pred_merged$preds.4 > mean(pred_merged$preds.4), "HR", "LR")
pred_merged$preds.bi.5 <- ifelse(pred_merged$preds.5 > mean(pred_merged$preds.5), "HR", "LR")
pred_merged$preds.bi.6 <- ifelse(pred_merged$preds.6 > mean(pred_merged$preds.6), "HR", "LR")
pred_merged$preds.bi.7 <- ifelse(pred_merged$preds.7 > mean(pred_merged$preds.7), "HR", "LR")
pred_merged$preds.bi.8 <- ifelse(pred_merged$preds.8 > mean(pred_merged$preds.8), "HR", "LR")
pred_merged$preds.bi.9 <- ifelse(pred_merged$preds.9 > mean(pred_merged$preds.9), "HR", "LR")
pred_merged$preds.bi.10 <- ifelse(pred_merged$preds.10 > mean(pred_merged$preds.10), "HR", "LR")
pred_merged$preds.bi.11 <- ifelse(pred_merged$preds.11 > mean(pred_merged$preds.11), "HR", "LR")
pred_merged$preds.bi.12 <- ifelse(pred_merged$preds.12 > mean(pred_merged$preds.12), "HR", "LR")
pred_merged$preds.bi.13 <- ifelse(pred_merged$preds.13 > mean(pred_merged$preds.13), "HR", "LR")
pred_merged$preds.bi.14 <- ifelse(pred_merged$preds.14 > mean(pred_merged$preds.14), "HR", "LR")
pred_merged$preds.bi.15 <- ifelse(pred_merged$preds.15 > mean(pred_merged$preds.15), "HR", "LR")

#  sampsInBothDatasets <- clinic_data[, "submitter_id"][clinic_data[, "submitter_id"] %in% newNames]
#  sampsInBothDatasets_TF <- clinic_data[, "submitter_id"] %in% newNames
#  for (j in 1:15){
#file  <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Surv/#Surv_",cate,"_",drug_names[i],"_",j,".pdf", sep="",collapse = NULL)
##title <- paste(cate," Survival Analysis *",drug_names[i]," * ",method_list[[j]], sep="",collapse = NULL)
#title <- method_list[[j]]
#var   <- paste("preds.bi",".",j, sep="",collapse = NULL)
#surv <- TCGAanalyze_survival(pred_merged,
#                     var,filename = file,
#                     xlab = title,height = 6, width=10,
#                   risk.table = FALSE,conf.int = F)
#}
#}
#library(pdftools)
#result_list <- list.files(path="/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Surv/",full.names #= TRUE)
#file  <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/#Surv_",cate,"_",drug_names[i],".pdf", sep="",collapse = NULL)
#result_list
#pdf_combine(result_list, output = file)

### ==================================================================
### == Survival analysing using TCGA data. Ref: https://github.com/BioinformaticsFMRP/TCGAbiolinks/blob/2c60481926592ce1a0c08c6229835d4eadd2b57c/R/methylation.R#L130
### ==================================================================
#library(survival)
#library(survminer)
#library(gridExtra)
#library(stringr)
#========================  Cleaning survival data  ========================
surv_plot <- list()
file  <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Surv/Surv_",cate,"_",drug_names[i],".pdf", sep="",collapse = NULL)
pdf(file,width=15, height=10)

for (k in 1:15){
name <- paste("preds.bi.",k, sep="",collapse = NULL)
pred_merged$pred <- as.factor(pred_merged[,name])
levels(pred_merged$pred) <- c("LR","HR")
pred_merged$gender <- as.factor(pred_merged$gender)
pred_merged$race <- as.factor(pred_merged$race)
data <- pred_merged[,c("submitter_id","drug_name","vital_status","days_to_death","days_to_last_follow_up","gender","age_at_index","race","pred")]

notDead <- is.na(data$days_to_death)

if (any(notDead == TRUE)) {
    data[notDead,"days_to_death"] <- data[notDead,"days_to_last_follow_up"]
}
# create a column to be used with survival package, info need
# to be TRUE(DEAD)/FALSE (ALIVE)
data$s <- grepl("dead|deceased",data$vital_status,ignore.case = TRUE)
#========================  Create formula for survival analysis  ========================
surv_model <- coxph(Surv(as.numeric(days_to_death),event=s)  ~ pred  + race + age_at_index, data)
#surv_model <- coxph(Surv(as.numeric(days_to_death),event=s)  ~ pred + gender + race + age_at_index, data)

surv_result0 <- data.frame(drug=drug_names[i],project=cate,N=nrow(data))
surv_result0$method <- method_list[k]
surv_result0$est.HR <- summary(surv_model)$coefficients[1,2]
surv_result0$p.value <- summary(surv_model)$coefficients[1,5]

surv_result <- rbind(surv_result,surv_result0)
#surv_result <- surv_result[surv_result$drug!="Gemcitabine",]
#surv_result$N <- ifelse(surv_result$drug=="Cisplatin",nrow(data),194)
#========================  Plotting for survival analysis  ========================
f.m <- formula(Surv(as.numeric(data$days_to_death),event=data$s) ~
                    data$pred)
fit <- do.call(survfit, list(formula = f.m, data = data))

#label.add.n <- function(x) {
#    paste0(x, " (n = ",
#           nrow(data[data[,"pred"] == x,]), ")")
#}
#
#
#    d <- survminer::surv_summary(fit, data = data)
##    order <- unname(sapply(levels(d$pred), function(x) unlist(str_split(x, "="))[2]))
#    order <- levels(d$pred)
#    labels <- sapply(order,label.add.n)
#
if(length(xlim) == 1) {
    xlim <- c(0,xlim)
}
    surv_plot[k] <-ggsurvplot(
        fit,                       # survfit object with calculated statistics.
        pval = paste("p=",signif(surv_result0$p.value,digits = 3),"\nestHR=",signif(surv_result0$est.HR,digits = 3),sep="",collapse = NULL),             # show p-value of log-rank test.
        xlab = method_list[[k]],               # customize X axis label.
        ggtheme = theme_light(),   # customize plot and risk table with a theme.
        palette =  rainbow(2)      # custom color palettes.
#        legend.labs = labels
    )

}

do.call("grid.arrange", c(surv_plot, ncol=5))
dev.off()

file  <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/Surv_stat/Surv_",cate,"_",drug_names[i],".csv", sep="",collapse = NULL)
write.csv(surv_result,file,quote = FALSE)
#suppressWarnings({
#    surv <- ggsurvplot(
#        fit,                       # survfit object with calculated statistics.
#        risk.table = risk.table,   # show risk table.
#        pval = pvalue,             # show p-value of log-rank test.
#        conf.int = conf.int,       # show confidence intervals for point estimaes of survival curves.
#        xlim = xlim,               # present narrower X axis, but not affect survival estimates.
#        main = main,               # Title
#        xlab = xlab,               # customize X axis label.
#        ggtheme = theme_light(),   # customize plot and risk table with a theme.
#        legend.title = legend,     # Legend title
#        legend.labs = labels,      # change legend labels.
#        palette =  color,          # custom color palettes.
#        ...
#    )
#})
}

}

source("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/2.3Validation-TCGA-Clinial.R")
