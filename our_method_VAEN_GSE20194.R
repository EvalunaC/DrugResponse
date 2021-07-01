source("Creating_drug_data.R")

load("/extraspace/ychen42/Drug_Response/VAEN/VAEN/Figure/Figure5/GSE20194/GSE20194.RData")
val.RPKM = expr.mat

#match("PLX-4720",possibleDrugs2)
drug <- "Paclitaxel"
drug_data <- getDrugData(drug)
trainFrame <- drug_data
