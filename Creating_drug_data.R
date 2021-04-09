set.seed(1000)

library(caret)

library(pRRophetic)
library(car)

library(glmnet)
library(dplyr)
library(magrittr)
library(MASS)
library(data.table)
library(parallel)
library(boot)
library(performance)
library(caret)
library(kernlab)
library(robustbase)
library(mgcv)

library(qrnn)



setwd("/extraspace/ychen42/Drug_Response/GROrignial/dataIn/rnaSeq/BRCA.Merge_rnaseqv2")
tpmDatMat_bc <- read.delim("BRCA.rnaseqv2_data.txt", as.is=T)

tpmDatMat_bc_tpm <- apply(tpmDatMat_bc[-1,which(tpmDatMat_bc[1,] == "scaled_estimate")], 2, as.numeric)
tpmDatMat_bc_tpm <- tpmDatMat_bc[-1,which(tpmDatMat_bc[1,] == "scaled_estimate")]
tpmDatMat_bc_tpm <- apply(tpmDatMat_bc_tpm, 2, as.numeric)
geneNames <- do.call(cbind, strsplit(tpmDatMat_bc[, "Hybridization.REF"], "|", fixed=TRUE))[1,][-1]
rownames(tpmDatMat_bc_tpm) <- geneNames
colnames(tpmDatMat_bc_tpm) <- substr(colnames(tpmDatMat_bc_tpm), 1, 28)
tpmDatMat_bc_tpm_logged <- log((tpmDatMat_bc_tpm*1000000)+1)

setwd("/extraspace/ychen42/Drug_Response/Data/")
CCLE_2018<-read.delim("CCLE_RNAseq_genes_rpkm_20180929.gct.txt") ##56202*1019
CCLE_2018<-CCLE_2018[,-1]

CCLE_2018_mat<-as.matrix(CCLE_2018)
rownames(CCLE_2018_mat)<-GeneNameCCLE<-CCLE_2018_mat[,1]
CCLE_2018_mat<-apply(CCLE_2018_mat[,-1],2,as.numeric)
rownames(CCLE_2018_mat)<-GeneNameCCLE
CCLE_2018_mat[1:3,1:3]
colnames(CCLE_2018_mat) <- gsub("^X", "",  colnames(CCLE_2018_mat))
cell_CCLE<- do.call(cbind, strsplit(colnames(CCLE_2018_mat), "_", fixed=TRUE))[1,] ## 1019
colnames(CCLE_2018_mat)<-cell_CCLE

GDSC2<-read.csv("GDSC2_IC50_matrix.csv", header=T,row.names=1)

#GDSC2<-read.csv("GDSC2_matrix.csv", header=T,row.names=1)
possibleDrugs2<-rownames(GDSC2) ##192
A2<-gsub("[[:punct:]]", " ", colnames(GDSC2))
colnames(GDSC2)<-gsub("[[:space:]]", "", A2) ##809

getCGPinfo_New<-function(drug){whichNas <- which(is.na(GDSC2[drug,]))
#drug_IC50<-GDSC_AUC_matrix[drug,][-whichNas]
drug_IC50<-GDSC2[drug,][-whichNas]


commonCellLines<-intersect(colnames(CCLE_2018_mat),colnames(drug_IC50))
CCLE_train<-CCLE_2018_mat[,commonCellLines]
drug_IC50_train<-drug_IC50[commonCellLines]

CCLE_train<-as.matrix(CCLE_train)
drug_IC50_train<-as.numeric(drug_IC50_train)
return(list(drug_IC50_train=drug_IC50_train, CCLE_train=CCLE_train))
}

######################################################
# Obtain medication list and prepare data from here
######################################################

getDrugData <- function(input_medication_name){

CCLETrainData<-getCGPinfo_New(input_medication_name)

#CCLETrainData<-getCGPinfo_New("Lapatinib")
#CCLETrainData<-getCGPinfo_New("Vinblastine")

#"Camptothecin"
#"Vinblastine"

trainingExprData<-CCLETrainData$CCLE_train
trainingPtype<-CCLETrainData$drug_IC50_train
testExprData<-tpmDatMat_bc_tpm_logged
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
return(trainFrame)  }



#input_medication_name <- "Vinblastine"
#input_medication_name <- "Camptothecin"
