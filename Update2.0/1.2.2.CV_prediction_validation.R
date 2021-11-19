#Read DL result on server1:
#========   Generate combined result
result_list <- list.files(path="/extraspace/ychen42/1.2.1.CV_perdiction_DL/")
result0 <- data.frame()
length(result_list)
  for (i in 1:length(result_list)){
    file <- paste("/extraspace/ychen42/1.2.1.CV_perdiction_DL/",result_list[i], sep="",collapse = NULL)
    result_file <- read.csv(file, header = FALSE,na.strings="NA")
    result0 <- rbind(result0,result_file)
  }
write.csv(result0, "/extraspace/ychen42/1.2.1.CV_perdiction_DL_combined.csv",quote = FALSE)

#Read other results on server2:
#========   Generate combined result

nonDL <- read.csv("/extraspace/ychen42/Drug_Response/Own_update2.0/Output/192CV_prediction/MK-2206.csv", header = TRUE,na.strings="NA")
DL <- read.csv("/extraspace/ychen42/Drug_Response/Own_update2.0/Output/1.2.1.CV_perdiction_DL_combined.csv", header = TRUE,na.strings="NA")
DL <- data.frame(value=DL$V1,y_test=DL$V2,drug_name=DL$V3,variable="Neural Network")

prediction <- data.frame()
result_list <- list.files(path="/extraspace/ychen42/Drug_Response/Own_update2.0/Output/192CV_prediction/")
length(result_list)
result_list <- strsplit(result_list, ".csv")

for (i in 1:length(result_list)){
  drug <- result_list[[i]]
file <- paste("/extraspace/ychen42/Drug_Response/Data/Drug192_TrainTest/",result_list[[i]],"_trainFrame.RData", sep="",collapse = NULL)
load(file)
drug_data <- trainFrame
  #========================================
  #=====    Split the data     ============
  #========================================
  train_length <- round(dim(drug_data)[1]*0.8)
  testResp <- drug_data[(train_length+1):dim(drug_data)[1],1]
  testResp_name <- rownames(drug_data)[(train_length+1):dim(drug_data)[1]]

   predict <- data.frame(y_test =testResp, id=testResp_name, drug_name=drug)
   prediction <- rbind(prediction,predict)
}

nonDL2 <- merge(nonDL, prediction, by=c("id","drug_name"),all.x=TRUE,all.y=FALSE, sort = TRUE)
nonDL <- nonDL2[c("value","y_test","drug_name","variable")]

data <- rbind(nonDL,DL)

data$R2 <- cor(data$value,data$y_test)
R2 <- data%>%
  group_by(drug_name,variable)%>%
  summarise(R2=cor(value,y_test))

stat <- data%>%
    group_by(drug_name,variable)%>%
    summarise(R2=cor(value,y_test), MAE=mean(abs(y_test - value)),y_test_mean=mean(y_test),value_mean=mean(value))

    data$y_cate <- ifelse(data$V1>=mean(data$V1),1,0)
    data$yhat_cate <- ifelse(data$V1.1>=mean(data$V1.1),1,0)
    (roc.pred2 <- roc(data$y_cate, data$yhat_cate, percent = TRUE, main = "Smoothing"))

write.csv(stat,"/extraspace/ychen42/Drug_Response/Own_update2.0/Output/stat.csv",quote = FALSE)

  library(ggplot2)

  pdf("/extraspace/ychen42/Drug_Response/Own_update2.0/Output/plots.pdf")
ggplot(stat, aes(x = R2, fill = variable)) +                       # Draw overlaying histogram
geom_histogram(position = "identity", alpha = 0.2, bins = 50)
  dev.off()
#library(dplyr)
#
#summary <- data%>%
#  group_by(drug_name,variable)%>%
#  summarise(value_mean=mean(value), y_test_mean=mean(y_test))
#data2 <- merge(data, summary, by=c("variable","drug_name"),all.x=TRUE,all.y=TRUE, sort = TRUE)
#
#data2 <- data %>%
#  group_by(drug_name, variable) %>%
#  mutate(bin = cut(value, breaks = c(-Inf, mean(value), Inf), labels = c(0,1)))
#
