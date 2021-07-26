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
