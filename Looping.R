



length(possibleDrugs2)
drug_data <- getDrugData(possibleDrugs2[1])
drug_result<- methods_result(drug_data, possibleDrugs2[1])
write.csv(df, "/extraspace/ychen42/Drug_Response/yiqings_work/Output/192Drug_MethodsResult.csv",
          col.names = FALSE,
          row.names = TRUE,
          quote = FALSE)
for (i in 2:length(possibleDrugs2)){
  drug_data <- getDrugData(possibleDrugs2[i])
  drug_result<- methods_result(drug_data, possibleDrugs2[i])
    write.csv(df, "/extraspace/ychen42/Drug_Response/yiqings_work/Output/192Drug_MethodsResult.csv",
          append = TRUE,
          col.names = TRUE,
          row.names = TRUE,
          quote = FALSE)
 
  }

#name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output/",drug,".csv", sep="",collapse = NULL)
#write.csv(df,name, quote=FALSE)




