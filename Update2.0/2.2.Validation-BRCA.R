#==========================================================================================================
#========   2.2 Validation using our training data and TCGA BRCA  ===========================
#==========================================================================================================
set.seed(1)

library(pRRophetic)
library(car)
library(ridge)



#source("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/0.1.Creating_drug_data.R")
#source("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/1.0.build_models_stateva.R")

load("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/Lapatinib_trainFrame.RData")
load("/extraspace/ychen42/Drug_Response/Data/BRCA.RData")



tpmDatMat_bc_tpm_logged <- BRCA

#Then the same as creating data.
CCLETrainData<-getCGPinfo_New("Lapatinib")
#> dim(CCLETrainData[[2]])
#[1] 56202   513
#> length(CCLETrainData[[1]])
#[1] 513
trainingExprData<-CCLETrainData$CCLE_train #[1] 56202   513
trainingPtype<-CCLETrainData$drug_IC50_train #513
testExprData<-tpmDatMat_bc_tpm_logged # is BRCA: [1] 20501  1104
selection=1
batchCorrect="standardize"
removeLowVaryingGenes=0.2
removeLowVaringGenesFrom="rawData"
printOutput=TRUE
powerTransformPhenotype=TRUE

trainExprMat=trainingExprData
testExprMat=testExprData


doVariableSelection <- function(exprMat, removeLowVaryingGenes)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}

homData <- homogenizeData(testExprData, trainingExprData,
                          batchCorrect = batchCorrect, selection = selection, printOutput = printOutput)
keepRows <- seq(1:nrow(homData$train))
evaluabeGenes <- rownames(homData$test)
keepRowsTrain <- doVariableSelection(trainingExprData[evaluabeGenes,
], removeLowVaryingGenes = removeLowVaryingGenes)
keepRowsTest <- doVariableSelection(testExprData[evaluabeGenes,
], removeLowVaryingGenes = removeLowVaryingGenes)
keepRows <- intersect(keepRowsTrain, keepRowsTest)
numberGenesRemoved <- nrow(homData$test) - length(keepRows)
if (printOutput)
  cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."))
offset = 0
if (powerTransformPhenotype) {
  if (min(trainingPtype) < 0) {
    offset <- -min(trainingPtype) + 1
    trainingPtype <- trainingPtype + offset
  }
  transForm <- powerTransform(trainingPtype)[[6]]
  trainingPtype <- trainingPtype^transForm
}

trainFrame <- data.frame(Resp = trainingPtype, t(homData$train[keepRows,]))
testFrame <- data.frame(t(homData$test[keepRows, ]))
#> dim(trainFrame)
#[1]  513 13226
#> dim(testFrame)
#[1]  1104 13225

#========================================
#  Then use  ==1.0.1. Building Models==
#========================================

preds <- preds_from_testFrame


##### Clinical outcome
# From "breast_cancer_analysis.html"
clinDataBrca <- read.delim("/extraspace/ychen42/Drug_Response/GROrignial/dataIn/clinical/nationwidechildrens.org_clinical_patient_brca.txt", as.is=T)
her2status <- clinDataBrca[, "her2_status_by_ihc"]
names(her2status) <- clinDataBrca[, "bcr_patient_barcode"]

sampleNames <- colnames(tpmDatMat_bc_tpm_logged)
theTumorSamples <- which(substring(sampleNames, 14, 16) == "01A") # identify the tumor samples, tumor samples annotated as "01" by TCGA, normal samples as "10".
newNames <- gsub(".", "-", substring(colnames(tpmDatMat_bc_tpm_logged), 1, 12), fixed=T)


pdf("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/BRCA_boxplot_clinical_15methods.pdf", width=8, height=10)
par(mfrow=c(3,5))
for (i in 1:15){
bcaPreds <- unlist(preds[[i]])

names(bcaPreds) <- newNames
bcaPreds <- bcaPreds[theTumorSamples] # Only include the tumor samples in this analysis. Results on normal samples are meaningless.
sampsInBothDatasets <- clinDataBrca[, "bcr_patient_barcode"][clinDataBrca[, "bcr_patient_barcode"] %in% newNames]
sampsInBothDatasets_TF <- clinDataBrca[, "bcr_patient_barcode"] %in% newNames

her2Neg <- which(her2status[sampsInBothDatasets_TF] == "Negative")
her2Pos <- which(her2status[sampsInBothDatasets_TF] == "Positive")
her2Equiv <- which(her2status[sampsInBothDatasets_TF] == "Equivocal")
(wiltest <- wilcox.test(bcaPreds[sampsInBothDatasets][her2Neg], bcaPreds[sampsInBothDatasets][her2Pos]))
#print(t.test(bcaPreds[sampsInBothDatasets][her2Neg], bcaPreds[sampsInBothDatasets][her2Pos]))

boxplot(list(Negative=bcaPreds[sampsInBothDatasets][her2Neg], Equivocal=bcaPreds[sampsInBothDatasets][her2Equiv], Positive=bcaPreds[sampsInBothDatasets][her2Pos]), las=1, col=c("#66c2a5", "#fc8d62", "#8da0cb"), pch=20, width=c(.75, .75, .75), ylab="Predicted Lapatinib Sensitivity", xlab=paste(method_list[i], "\n Wilcxn RS p = ", signif(wiltest$p.value,digits = 3), sep=""),cex.axis=.75, outcol="#00000033")
}
dev.off()
