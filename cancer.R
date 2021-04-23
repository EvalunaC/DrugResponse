#' I now want to investigate associations between mutations and between drug sensitivty when I look across the entire TCGA datasets. I will summarize the mutations on a per-gene basis, i.e. I will call a gene "mutated" if there is any protein coding change anywhere in the gene.

#' Load the somatic mutation data from TCGA, this data is from exome sequencing and was obtained from firebrowse.org.
#' First I need to create a matrix of genes by samples and indicate whether samples do / don't have somatic mutations. I.e. summarize the mutation data by gene.
theRootDir <- "/mnt/data_scratch/finalData/"
dir.create(paste(theRootDir, "tables/", sep=""), showWarnings = FALSE)
dirList <- dir(paste(theRootDir, "dataIn/mutation_data/", sep=""))
subfolders <- dirList[-grep(".tar.gz", dirList)]
mutsListAll <- list() # a list of mutations occuring in each sample.
allFileNames <- character()
proteinChangingMutations <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "In_Frame_Del", "Frame_Shift_Ins", "In_Frame_Ins", "Nonstop_Mutation", "De_novo_Start_OutOfFrame", "De_novo_Start_InFrame", "Missense", "Read-through", "Indel")
genesWithProteinChangeList <- list()
for(i in 1:length(subfolders))
{
  subFolderFiles <- dir(paste(theRootDir, "dataIn/mutation_data/", subfolders[i], sep=""))[-1] # list the files (getting rid of manifest.txt)
  tcgaIds <- sapply(subFolderFiles, function(item)return(strsplit(item, ".", fixed=T)[[1]][1]))
  allFileNames <- c(allFileNames, subFolderFiles)
  mutsListAll[[i]] <- list()
  genesWithProteinChangeList[[i]] <- list()

  # now for each file in this subfolder, pull out a list of genes with somatic mutations
  for(j in 1:length(subFolderFiles))
  {
    mutfile <- read.delim(paste(theRootDir, "dataIn/mutation_data/", subfolders[i], "/", subFolderFiles[j], sep=""), as.is=T)
    variantType <- mutfile[ ,"Variant_Classification"]
    theGenes <- mutfile[ ,"Hugo_Symbol"]
    names(variantType) <- theGenes
    mutsListAll[[i]][[j]] <- variantType
    genesWithProteinChangeList[[i]][[j]] <- unique(names(variantType[variantType %in% proteinChangingMutations]))
  }
}
allMutatedGenes <- unique(names(unlist(mutsListAll)))
mutationTypes <- table(unlist(mutsListAll))
mutsListAll_unlist <- unlist(mutsListAll, recursive=F)
genesWithProteinChangeList_unlist <- unlist(genesWithProteinChangeList, recursive=F)

#' List the frequency of the different types of mutations catalogued in the files above.
print(sort(mutationTypes, decreasing=T))


#' From mutsListAll we can then create a matrix indicating if the gene has a coding mutation.
tcgaIdsAll <- sapply(strsplit(allFileNames, ".", fixed=T), function(l)return(l[[1]]))
mutMat <- numeric((length(tcgaIdsAll)*length(allMutatedGenes)))
dim(mutMat) <- c(length(allMutatedGenes), length(tcgaIdsAll))
rownames(mutMat) <- allMutatedGenes
colnames(mutMat) <- tcgaIdsAll
print(mutMat[1:5, 1:5])
mutMatTest <- mutMat

#' Now populate this matrix with the relevant information about what kind of mutatation each gene has in each sample.
for(i in 1:length(tcgaIdsAll))
{
  mutMat[genesWithProteinChangeList_unlist[[i]], i] <- rep(1, length(genesWithProteinChangeList_unlist[[i]]))
}
print(mutMat[1:5, 1:5])


#' What are the most commonly mutated genes, is this consistent with expectation? Yes.
numMuts <- apply(mutMat, 1, function(row)return(sum(!row == 0)))
print(sort(numMuts, decreasing=T)[1:10])


#' These TCGA samples contain several different types of tumors. Here, we will focus on the largest of these groups ("Primary Solid Tumors"). These are 8528 samples and are labeled "01" in the TCGA sample IDs.
#' 01: Primary Solid Tumor
#' 02: Recurrent Solid Tumor
#' 03: Primary Blood derived cancer - peripheral blood
#' 04: Regcurrent blood rerived cancer - bone marrow
#' 05: Additional - new primary
#' 06: Metastatic
tumorTypeId <- sapply(strsplit(colnames(mutMat), "-", fixed=TRUE), function(l)return(l[4]))
print(table(tumorTypeId))

#' Lets remove everything but the "Primary Solid Tumors (i.e. "01")".
mutMat_only01 <- mutMat[, tumorTypeId == "01"]
theIds <- colnames(mutMat_only01)
mutMat_nodups <- mutMat_only01[, !duplicated(theIds)] # Some samples were listed in multiple folders when the data was downloaded from TCGA (e.g. the KIRC and KIPAN folders). This step will remove these samples.
mutIds <- sapply(strsplit(colnames(mutMat_nodups), "-", fixed=T), function(l)return(l[3]))
colnames(mutMat_nodups) <- mutIds

# save(mutMat_nodups, file="/mnt/data_scratch/prediXcanProj/Results/rDatas/mutMat.RData")

#' Load the imputed drug sensitivty data.
load(file=paste(theRootDir, "dataOut/allDrugPredictions_mat.RData", sep="")) # allDrugPredictions_mat, cancerTypesVec,
names(cancerTypesVec) <- colnames(allDrugPredictions_mat)

#' Make a plot for lapatinib across all cancer types
lapatinibTypes <- split(allDrugPredictions_mat["Lapatinib", ], cancerTypesVec)
svg(paste(theRootDir, "figures/lapatinibTypes.svg", sep=""), width=8, height=8)
boxplot(lapatinibTypes, las=2)
dev.off()

#' Extract the 01a samples, i.e. tumor samples.
all01ASamples <- colnames(allDrugPredictions_mat)[which(sapply(strsplit(colnames(allDrugPredictions_mat), ".", fixed=T), function(a)a[4]) == "01A")]
preds01a <- allDrugPredictions_mat[, all01ASamples]
cancerTypes01a <- cancerTypesVec[all01ASamples]
sampIds01a <- sapply(strsplit(all01ASamples, ".", fixed=T), function(l)return(l[3]))
names(cancerTypes01a) <- sampIds01a
colnames(preds01a) <- sampIds01a
inPredAndMutData <- sampIds01a[sampIds01a %in% mutIds] # samples for which we have both predicted drug response and mutation calls

#' Run the associations between all genes and drugs, for drugs with at least 50 mutations.
cancerTypes01a_filt_ord <- cancerTypes01a[inPredAndMutData]
preds01a_filt_ord <- preds01a[, inPredAndMutData]
mutMat_nodups_ordFilt <- mutMat_nodups[, inPredAndMutData]
commonMuts <- apply(mutMat_nodups_ordFilt, 1, sum)
commonlyMutated <- mutMat_nodups_ordFilt[which(commonMuts >= 50), ]
commonlyMutated <- commonlyMutated[-which(rownames(commonlyMutated) == "Unknown"), ] # there is an entry for gene wiht an "Unknown" HUGO ID. Remove these.
pValList <- list()
betaValList <- list()
for(i in 1:nrow(preds01a_filt_ord))
{
  pValList[[i]] <- numeric()
  betaValList[[i]] <- numeric()
  for(j in 1:nrow(commonlyMutated))
  {
    thecoefs <- coef(summary(lm(preds01a_filt_ord[i,]~commonlyMutated[j,]+cancerTypes01a_filt_ord)))
    pValList[[i]][[j]] <- thecoefs[2,4]
    betaValList[[i]][[j]] <- thecoefs[2,1]
  }
}


#' Create the "cancer-type" corrected TCGA drug prediction data matrix. This can potentially be used by others in subsequent analysis. This will be Supplementary Table 9
imputedDrugResponses_correctCantype <- numeric(nrow(preds01a_filt_ord)*ncol(preds01a_filt_ord))
dim(imputedDrugResponses_correctCantype) <- c(nrow(preds01a_filt_ord), ncol(preds01a_filt_ord))
for(i in 1:nrow(preds01a_filt_ord))
{
  imputedDrugResponses_correctCantype[i, ] <- residuals(lm(preds01a_filt_ord[i,]~cancerTypes01a_filt_ord)) ########为啥取的residuals？ cancerTypes01a_filt_ord 是cancer type。
}
rownames(imputedDrugResponses_correctCantype) <- rownames(preds01a_filt_ord)
colnames(imputedDrugResponses_correctCantype) <- colnames(preds01a_filt_ord)
write.csv(imputedDrugResponses_correctCantype, paste(theRootDir, "tables/imputedDrugResponses_correctCantype.csv", sep="")) #table 9


#' The total number of genes genome wide with at least one somatic protein coding change called in TCGA.
length(commonMuts) # [1] 27856

#' The number of genes which have a somatic protein coding change in at least 50 samples
print(sum(commonMuts > 50)) # [1] 1673

# Get the adjusted p-value for each gene-drug combination, pull out the significant associations and create a supplementary table that lists these for "predictable" drugs?.
sigPs <- list()
pAdjListCantype <- list()
for(i in 1:length(pValList))
{
  names(pValList[[i]]) <- rownames(commonlyMutated)
  names(betaValList[[i]]) <- rownames(commonlyMutated)
  padj <- p.adjust(pValList[[i]], method="BH")
  sigPs[[i]] <- padj[padj < 0.05]
  pAdjListCantype[[i]] <- padj
}
names(sigPs) <- rownames(preds01a_filt_ord)
names(pValList) <- rownames(preds01a_filt_ord)
names(betaValList) <- rownames(preds01a_filt_ord)
names(pAdjListCantype) <- rownames(preds01a_filt_ord)

#' Print the top associations
print(sort(unlist(pValList))[1:30])

#' Write all significant associations out to a supplmentary table. We are only considering results for "predictable" drugs, i.e. drugs with a spearman corrleation of > 0.3 in cross validation in the GDSC cancer cell lines on which these models were intially fit.
predictableDrugs <- c("ABT.263", "ABT.888", "AG.014699", "AICAR", "ATRA", "Axitinib", "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", "AZD7762", "BAY.61.3606", "BIBW2992", "Bicalutamide", "BI.D1870", "Bleomycin", "BMS.536924", "BMS.754807", "Bortezomib", "Bosutinib", "BX.795", "Camptothecin", "CEP.701", "CGP.082996", "CHIR.99021", "CI.1040", "Cisplatin", "Cytarabine", "Dasatinib", "DMOG", "Docetaxel", "Elesclomol", "Erlotinib", "FH535", "FTI.277", "Gefitinib", "Gemcitabine", "IPA.3", "Lapatinib", "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", "Nilotinib", "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", "PAC.1", "PD.0325901", "PD.0332991", "PD.173074", "PLX4720", "RDEA119", "SB590885", "Sunitinib", "Temsirolimus", "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", "Vorinostat", "VX.702", "WH.4.023", "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", "ZM.447439")
print(length(predictableDrugs))

# We only want to report the p-values for the "predictable" drugs....
pValList_predictable <- unlist(pValList[predictableDrugs])
ord <- order(pValList_predictable)
pValList_predictable_ord <- pValList_predictable[ord]
betaValList_predictable <- unlist(betaValList[predictableDrugs])[ord]
pAdjListCantype_predictable <- unlist(pAdjListCantype[predictableDrugs])[ord]
pAdjListCantype_predictable_forNumModels <- (pAdjListCantype_predictable*71)
pAdjListCantype_predictable_forNumModels[pAdjListCantype_predictable_forNumModels > 1] <- 1 # adjust the FDRs for the number of models and cap this at 1. This is equivalent to a Bonferroni correction for 71 tests.
outTab <- cbind(pValList_predictable_ord, pAdjListCantype_predictable, betaValList_predictable, pAdjListCantype_predictable_forNumModels)
colnames(outTab) <- c("P-value", "FDR", "Effect Size (beta)", "FDR Corrected for 71 Models")
write.csv(outTab, paste(theRootDir, "tables/allResults_control_for_tissue.csv", sep=""))
sum(outTab[,4] < 0.05) # [1] 104

#' Make figures for Nultlin 3 and P53, for PD.0332991 and RB1 and for KRAS and Erlotinib.
svg(paste(theRootDir, "figures/nutlin3_p53.svg", sep=""), width=3, height=3)
hist(-log10(pValList[["Nutlin.3a"]]), breaks=100, col="#8dd3c7", las=1, xlab=expression("-Log"[10]*"P-value"), cex.axis=0.75, main="Nutlin-3a")
abline(v=-log10(pValList[["Nutlin.3a"]]["TP53"]), col="red")
dev.off()

svg(paste(theRootDir, "figures/PD0332991_rb1.svg", sep=""), width=3, height=3)
hist(-log10(pValList[["PD.0332991"]]), breaks=100, col="#8dd3c7", main="PD-0332991", las=1, xlab=expression("-Log"[10]*"P-value"), cex.axis=0.75)
abline(v=-log10(pValList[["PD.0332991"]]["RB1"]), col="red")
dev.off()

svg(paste(theRootDir, "figures/erlotinib_kras.svg", sep=""), width=3, height=3)
hist(-log10(pValList[["Erlotinib"]]), breaks=100, col="#8dd3c7", main="Erlotinib", las=1, xlab=expression("-Log"[10]*"P-value"), cex.axis=0.75)
abline(v=-log10(pValList[["Erlotinib"]]["KRAS"]), col="red")
dev.off()



#' Now re-do the gene-drug association analysis, but control for "general levels of drug sensitivty (glds)" instead of for cancer type. This algorithm is described in Geeleher et al, Genome Biolgoy (2016).
#' First, I need to calculate a set of negative control (i.e. unrelated) drugs for each drug. Then I need to calculate the principal components on the IC50 values for these drugs.
drugRelatedness <- read.csv(paste(theRootDir, "dataIn/categorys.csv", sep=""), as.is=TRUE)
theseDrugNames <- rownames(preds01a_filt_ord)
drugRelatedness[, "theseDrugNames"] <- unlist(strsplit(drugRelatedness[, "DrugNamesOtherFile"], "_IC_50"))
pairCor <- cor(t(preds01a_filt_ord), method="spearman")#####每两个药取correlation
controlPcsList <- list()
for(j in 1:nrow(preds01a_filt_ord))
{
  categoryThisDrug <- drugRelatedness[, "Drug.Category"][which(drugRelatedness["theseDrugNames"] == rownames(preds01a_filt_ord)[j])]
  negControlDrugs <- na.omit(drugRelatedness[!drugRelatedness[, "Drug.Category"] %in% categoryThisDrug, "theseDrugNames"])
  pairwiseCorNear <- names(rank(abs(pairCor[, colnames(t(preds01a_filt_ord))[j]]))[118:137]) # also remove very correlated drugs...
  negControlDrugs <- setdiff(negControlDrugs, pairwiseCorNear)
  controlPCsAll <- prcomp(t(preds01a_filt_ord)[, negControlDrugs])$x
  controlPcsList[[j]] <- controlPCsAll
}

######################################################################################################################################################################

#' Run the analysis and control for GLDS, for each drug, control for the first 50 principal components of the negative control drugs.
pValList_glds_only <- list()
betaValList_glds_only <- list()
for(i in 1:nrow(preds01a_filt_ord))
{
  pValList_glds_only[[i]] <- numeric()
  betaValList_glds_only[[i]] <- numeric()
  for(j in 1:nrow(commonlyMutated))
  {
    thecoefs <- coef(summary(lm(preds01a_filt_ord[i,]~commonlyMutated[j,]+controlPcsList[[i]][,1:50])))
    pValList_glds_only[[i]][[j]] <- thecoefs[2,4]
    betaValList_glds_only[[i]][[j]] <- thecoefs[2,1]
  }
}

#' Create the "cancer-type" corrected TCGA drug prediction data matrix. This can potentially be used by others in subsequent analysis. This will be Supplementary Table 10
imputedDrugResponses_correctGlds <- numeric(nrow(preds01a_filt_ord)*ncol(preds01a_filt_ord))
dim(imputedDrugResponses_correctGlds) <- c(nrow(preds01a_filt_ord), ncol(preds01a_filt_ord))
for(i in 1:nrow(preds01a_filt_ord))
{
  imputedDrugResponses_correctGlds[i, ] <- residuals(lm(preds01a_filt_ord[i,]~controlPcsList[[i]][,1:50]))
}
rownames(imputedDrugResponses_correctGlds) <- rownames(preds01a_filt_ord)
colnames(imputedDrugResponses_correctGlds) <- colnames(preds01a_filt_ord)
write.csv(imputedDrugResponses_correctGlds, paste(theRootDir, "tables/imputedDrugResponses_correctGlds.csv", sep=""))


#' Calculate false discovery rates for these samples. Then output a table of the top results, but only include the results for the "predictable" drugs, i.e. those with a spearman's correlation of > 0.3 in the cell lines on which these data were fit.
sigPs_glds_only <- list()
pAdjListGlds <- list()
for(i in 1:length(pValList_glds_only))
{
  names(pValList_glds_only[[i]]) <- rownames(commonlyMutated)
  names(betaValList_glds_only[[i]]) <- rownames(commonlyMutated)
  padj <- p.adjust(pValList_glds_only[[i]], method="BH")
  sigPs_glds_only[[i]] <- padj[padj < 0.0005]
  pAdjListGlds[[i]] <- padj
}
names(sigPs_glds_only) <- rownames(preds01a_filt_ord)
names(pValList_glds_only) <- rownames(preds01a_filt_ord)
names(betaValList_glds_only) <- rownames(preds01a_filt_ord)
names(pAdjListGlds) <- rownames(preds01a_filt_ord)
print(sort(unlist(pValList_glds_only))[1:30])

#' We want to only cause adjusted p-values for the "predictable" drugs...
pValList_glds_only_predictable <- unlist(pValList_glds_only[predictableDrugs])
ord <- order(pValList_glds_only_predictable)
pValList_glds_only_predictable_ord <- pValList_glds_only_predictable[ord]
betaValList_glds_only_predictable <- unlist(betaValList_glds_only[predictableDrugs])[ord]
pAdjListGlds_predictable <- unlist(pAdjListGlds[predictableDrugs])[ord]
pAdjListGlds_predictable_forNumModels <- (pAdjListGlds_predictable*71)
pAdjListGlds_predictable_forNumModels[pAdjListGlds_predictable_forNumModels > 1] <- 1 # adjust the FDRs for the number of models and cap this at 1. This is equivalent to a Bonferroni correction for 71 tests.
outTab_glds <- cbind(pValList_glds_only_predictable_ord, pAdjListGlds_predictable, betaValList_glds_only_predictable, pAdjListGlds_predictable_forNumModels)
colnames(outTab_glds) <- c("P-value", "FDR", "Effect Size (beta)", "FDR corrected for number of models")
write.csv(outTab_glds, paste(theRootDir, "tables/allResults_glds.csv", sep=""))
sum(outTab_glds[,4] < 0.05) # [1] 263

#' Make soem figures for the top results.
gldspValMat <- do.call(cbind, pValList_glds_only)

#' Make GLDS figures for Nultlin 3 and P53, for PD.0332991 and RB1 and for KRAS and Erlotinib.
svg(paste(theRootDir, "figures/nutlin3_p53_glds.svg", sep=""), width=3, height=3)
hist(-log10(pValList_glds_only[["Nutlin.3a"]]), breaks=100, col="#8dd3c7", main="Nutlin-3a", las=1, xlab=expression("-Log"[10]*"P-value"), cex.axis=0.75)
abline(v=-log10(pValList_glds_only[["Nutlin.3a"]]["TP53"]), col="red")
dev.off()

svg(paste(theRootDir, "figures/PD0332991_rb1_glds.svg", sep=""), width=3, height=3)
hist(-log10(pValList_glds_only[["PD.0332991"]]), breaks=100, col="#8dd3c7", main="PD-0332991", las=1, xlab=expression("-Log"[10]*"P-value"), cex.axis=0.75)
abline(v=-log10(pValList_glds_only[["PD.0332991"]]["RB1"]), col="red")
dev.off()

svg(paste(theRootDir, "figures/erlotinib_kras_glds.svg", sep=""), width=3, height=3)
hist(-log10(pValList_glds_only[["Erlotinib"]]), breaks=100, col="#8dd3c7", main="Erlotinib", las=1, xlab=expression("-Log"[10]*"P-value"), cex.axis=0.75)
abline(v=-log10(pValList_glds_only[["Erlotinib"]]["KRAS"]), col="red")
dev.off()

#' Print some of the results for mutations that we'd expect for drugs which are "predictable". And for drugs which are not....
print(pValList_glds_only[["PD.0332991"]]["RB1"])
print(pValList_glds_only[["Nutlin.3a"]]["TP53"])
print(pValList_glds_only[["PLX4720"]]["BRAF"])
print(pValList_glds_only[["SB590885"]]["BRAF"])
print(pValList_glds_only[["PD.0325901"]]["BRAF"])
print(pValList_glds_only[["Gefitinib"]]["EGFR"])
print(pValList_glds_only[["Erlotinib"]]["KRAS"])
print(pValList_glds_only[["PD.0332991"]]["CDKN2A"])
print(pValList_glds_only[["AZD.2281"]]["BRCA2"]) # There are potential problems with screening of PARP inhibors in GDSC, see Heitmann et al (Oral Oncology 2014). Drugs likely not screened for long enough.
print(pValList_glds_only[["AZD.2281"]]["BRCA1"])
print(pValList_glds_only[["ATRA"]]["RARA"]) # not enough mutations
print(pValList_glds_only[["Sunitinib"]]["PDGFRB"])
print(pValList_glds_only[["Sunitinib"]]["FLT3"]) # not enough mutations
print(pValList_glds_only[["Erlotinib"]]["EGFR"])
print(pValList_glds_only[["Rapamycin"]]["AKT1"]) # Rapamycin not in the "predictable" drugs list
print(pValList_glds_only[["Imatinib"]]["KIT"])
print(pValList_glds_only[["Nilotinib"]]["KIT"])
print(pValList_glds_only[["Sunitinib"]]["KIT"])

#' Print the same associations when controlling for cancer type instead of GLDS. Arguably the GLDS results are slightly better.
print(pValList[["PD.0332991"]]["RB1"])
print(pValList[["Nutlin.3a"]]["TP53"])
print(pValList[["PLX4720"]]["BRAF"])
print(pValList[["SB590885"]]["BRAF"])
print(pValList[["PD.0325901"]]["BRAF"])
print(pValList[["Gefitinib"]]["EGFR"])
print(pValList[["Erlotinib"]]["KRAS"])
print(pValList[["PD.0332991"]]["CDKN2A"])
print(pValList[["AZD.2281"]]["BRCA2"])
print(pValList[["AZD.2281"]]["BRCA1"])
print(pValList[["ATRA"]]["RARA"]) # not enough mutations
print(pValList[["Sunitinib"]]["PDGFRB"])
print(pValList[["Sunitinib"]]["FLT3"]) # not enough mutations
print(pValList[["Erlotinib"]]["EGFR"])
print(pValList[["Rapamycin"]]["AKT1"]) # Rapamycin not in the "predictable" drugs list
print(pValList[["Imatinib"]]["KIT"])
print(pValList[["Sunitinib"]]["KIT"])
print(pValList[["Nilotinib"]]["KIT"])
print(pValList_glds_only[["Sunitinib"]]["KIT"])
