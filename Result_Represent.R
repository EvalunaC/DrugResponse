result_list <- list.files(path="/extraspace/ychen42/Drug_Response/yiqings_work/Output/")

result0 <- data.frame()

length(result_list)
for (i in 1:length(result_list)){
  file <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output/",result_list[i], sep="",collapse = NULL)
  result_file <- read.csv(file, header = TRUE,na.strings="NA")
  result0 <- rbind(result0,result_file)
}
#name <- paste("/extraspace/ychen42/Drug_Response/yiqings_work/Output/",drug,".csv", sep="",collapse = NULL)
#write.csv(df,name, quote=FALSE)
result_Cisplatin <- read.csv("/extraspace/ychen42/Drug_Response/yiqings_work/Output/Cisplatin.csv", header = TRUE,na.strings="NA")
result_Gefitinib <- read.csv("/extraspace/ychen42/Drug_Response/yiqings_work/Output/Gefitinib.csv", header = TRUE,na.strings="NA")
result <- rbind(result_Cisplatin,result_Gefitinib)

method_1 <- result[result$X==1,]

write.csv(result0, "/extraspace/ychen42/Drug_Response/yiqings_work/57Drug_MethodsResult.csv",          col.names = TRUE,
          row.names = TRUE,
          quote = FALSE)

mean(result$MAE, na.rm = TRUE)

##############################################
# Method performance (MAE, RMSE, R square)
##############################################
library(ggplot2)
pdf("/extraspace/ychen42/Drug_Response/yiqings_work/method_perform_violin.pdf")
p <- ggplot(result0, aes(factor(method), RMSE))
p <- p + geom_violin(aes(colour = "#1268FF"), alpha = 0.3)
q <- p + geom_violin(aes(y = R_Squared, colour = "#3268FF"), alpha = 0.3)
dev.off()

