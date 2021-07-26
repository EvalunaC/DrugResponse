setwd("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/")

source("0.1.Creating_drug_data.R")
source("/extraspace/ychen42/Drug_Response/Own_update2.0/some_code/1.0.build_models_stateva.R")


#==================================================
#========   1.1.1. Loop build_models by drug (N=192)
#==================================================

library(doMC)
options(cores = 25)
registerDoMC()

#drug="Lapatinib"
drug="NVP-ADW742"
drug_data <- getDrugData(drug)

169
for (i in 1:length(possibleDrugs2)){
  drug=possibleDrugs2[i]
  drug_data <- getDrugData(drug)[[1]]
  build_stateval(drug=drug,trainFrame=drug_data,testFrame=drug_data[,-1],outFolder="/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/4Stats_eachdrug/")
}

for (i in 169:length(possibleDrugs2)){
  drug=possibleDrugs2[i]
  drug_data <- getDrugData(drug)
  build_stateval(drug=drug,trainFrame=drug_data[[1]],testFrame=drug_data[[2]])
}


#==================================================
#========   1.1.2. Generate combined result
#==================================================
result_list <- list.files(path="/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/4Stats_eachdrug/")
result0 <- data.frame()
length(result_list)
  for (i in 1:length(result_list)){
    file <- paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/4Stats_eachdrug/",result_list[i], sep="",collapse = NULL)
    result_file <- read.csv(file, header = TRUE,na.strings="NA")
    result0 <- rbind(result0,result_file)
  }
write.csv(result0, "/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/192drugs_4Stats.csv",quote = FALSE)
