#source("rework.R") get everything and continue from rework

load("/extraspace/ychen42/Drug_Response/Data/brcaTestFrame.RData")

### For TCGA BRCA, only Lapatinib.
set.seed(1)

drug_data <- getDrugData("Lapatinib")
trainFrame <- drug_data

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


cv_brca <- c()
preds <- list()
for (i in 1:7){
  cat(i)
  preds[[i]]<-predict(methods_model[i],testFrame)
}
preds[[8]] <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame))
preds[[9]]  <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame))



##### Clinical outcome
# From "breast_cancer_analysis.html"
clinDataBrca <- read.delim("/extraspace/ychen42/Drug_Response/GROrignial/dataIn/clinical/nationwidechildrens.org_clinical_patient_brca.txt", as.is=T)
her2status <- clinDataBrca[, "her2_status_by_ihc"]
names(her2status) <- clinDataBrca[, "bcr_patient_barcode"]

sampleNames <- colnames(tpmDatMat_bc_tpm_logged)
theTumorSamples <- which(substring(sampleNames, 14, 16) == "01A") # identify the tumor samples, tumor samples annotated as "01" by TCGA, normal samples as "10".
newNames <- gsub(".", "-", substring(colnames(tpmDatMat_bc_tpm_logged), 1, 12), fixed=T)


pdf("/extraspace/ychen42/Drug_Response/yiqings_work/BRCA_boxplot_clinical_9methods.pdf", width=8, height=10)
par(mfrow=c(3, 3))
for (i in 1:9){
bcaPreds <- unlist(preds[[i]])

names(bcaPreds) <- newNames
bcaPreds <- bcaPreds[theTumorSamples] # Only include the tumor samples in this analysis. Results on normal samples are meaningless.
sampsInBothDatasets <- clinDataBrca[, "bcr_patient_barcode"][clinDataBrca[, "bcr_patient_barcode"] %in% newNames]
sampsInBothDatasets_TF <- clinDataBrca[, "bcr_patient_barcode"] %in% newNames

her2Neg <- which(her2status[sampsInBothDatasets_TF] == "Negative")
her2Pos <- which(her2status[sampsInBothDatasets_TF] == "Positive")
her2Equiv <- which(her2status[sampsInBothDatasets_TF] == "Equivocal")
(wiltest <- wilcox.test(bcaPreds[sampsInBothDatasets][her2Neg], bcaPreds[sampsInBothDatasets][her2Pos]))
print(t.test(bcaPreds[sampsInBothDatasets][her2Neg], bcaPreds[sampsInBothDatasets][her2Pos]))

boxplot(list(Negative=bcaPreds[sampsInBothDatasets][her2Neg], Equivocal=bcaPreds[sampsInBothDatasets][her2Equiv], Positive=bcaPreds[sampsInBothDatasets][her2Pos]), las=1, col=c("#66c2a5", "#fc8d62", "#8da0cb"), pch=20, width=c(.75, .75, .75), ylab="Predicted Lapatinib Sensitivity", xlab=paste(methods[i], "\n Wilcoxon rank sum p = ", signif(wiltest$p.value,digits = 3), sep=""),cex.axis=.75, outcol="#00000033")
}
dev.off()



pdf("/extraspace/ychen42/Drug_Response/yiqings_work/BRCA_boxplot_clinical_ranger.pdf")
seeds<- list(c(1,1,1,1,1,1), c(1,1,1,1,1,1))
model_ranger<-train(Resp~.,data=trainFrame,num.trees = 50,method="ranger",trControl = trainControl(method = "oob", seed = seeds))

bcaPreds<-predict(model_ranger,testFrame)

bcaPreds <- bcaPreds[theTumorSamples] # Only include the tumor samples in this analysis. Results on normal samples are meaningless.
sampsInBothDatasets <- clinDataBrca[, "bcr_patient_barcode"][clinDataBrca[, "bcr_patient_barcode"] %in% newNames]
sampsInBothDatasets_TF <- clinDataBrca[, "bcr_patient_barcode"] %in% newNames

her2Neg <- which(her2status[sampsInBothDatasets_TF] == "Negative")
her2Pos <- which(her2status[sampsInBothDatasets_TF] == "Positive")
her2Equiv <- which(her2status[sampsInBothDatasets_TF] == "Equivocal")
(wiltest <- wilcox.test(bcaPreds[sampsInBothDatasets][her2Neg], bcaPreds[sampsInBothDatasets][her2Pos]))
print(t.test(bcaPreds[sampsInBothDatasets][her2Neg], bcaPreds[sampsInBothDatasets][her2Pos]))

boxplot(list(Negative=bcaPreds[sampsInBothDatasets][her2Neg], Equivocal=bcaPreds[sampsInBothDatasets][her2Equiv], Positive=bcaPreds[sampsInBothDatasets][her2Pos]), las=1, col=c("#66c2a5", "#fc8d62", "#8da0cb"), pch=20, width=c(.75, .75, .75), ylab="Predicted Lapatinib Sensitivity", xlab=paste(methods[i], "\n Wilcoxon rank sum p = ", signif(wiltest$p.value,digits = 3), sep=""),cex.axis=.75, outcol="#00000033")

dev.off()
