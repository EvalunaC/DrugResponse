pred <- read.csv("/extraspace/ychen42/Drug_Response/DeepLearning/ForgeNet_predict.csv",header=FALSE)

##### Clinical outcome
# From "breast_cancer_analysis.html"
clinDataBrca <- read.delim("/extraspace/ychen42/Drug_Response/GROrignial/dataIn/clinical/nationwidechildrens.org_clinical_patient_brca.txt", as.is=T)
her2status <- clinDataBrca[, "her2_status_by_ihc"]
names(her2status) <- clinDataBrca[, "bcr_patient_barcode"]

sampleNames <- colnames(tpmDatMat_bc_tpm_logged)
theTumorSamples <- which(substring(sampleNames, 14, 16) == "01A") # identify the tumor samples, tumor samples annotated as "01" by TCGA, normal samples as "10".
newNames <- gsub(".", "-", substring(colnames(tpmDatMat_bc_tpm_logged), 1, 12), fixed=T)



bcaPreds <- pred

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
