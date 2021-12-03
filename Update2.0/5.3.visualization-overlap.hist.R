library(ggplot2)
library(reshape2)
library(dplyr)

PCC_Combine_CTRP1 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/PCC_Combine_CTRP1.csv",header=TRUE)
PCC_Combine_CTRP2 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/PCC_Combine_CTRP2.csv",header=TRUE)
PCC_Combine_GDSC1 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/PCC_Combine_GDSC1.csv",header=TRUE)
PCC_Combine_GDSC2 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/PCC_Combine_GDSC2.csv",header=TRUE)

PCC_Combine_CTRP1 <- data.frame(method=PCC_Combine_CTRP1$method,drug=PCC_Combine_CTRP1$drug,variable="PCC",value=PCC_Combine_CTRP1$PCC_CV28,dataset="CTRP1")
PCC_Combine_CTRP2 <- data.frame(method=PCC_Combine_CTRP2$method,drug=PCC_Combine_CTRP2$drug,variable="PCC",value=PCC_Combine_CTRP2$PCC_CV28,dataset="CTRP2")
PCC_Combine_GDSC1 <- data.frame(method=PCC_Combine_GDSC1$method,drug=PCC_Combine_GDSC1$drug,variable="PCC",value=PCC_Combine_GDSC1$PCC_CV28,dataset="GDSC1")
PCC_Combine_GDSC2 <- data.frame(method=PCC_Combine_GDSC2$method,drug=PCC_Combine_GDSC2$drug,variable="PCC",value=PCC_Combine_GDSC2$PCC_CV28,dataset="GDSC2")

pcc <- rbind(PCC_Combine_CTRP1,PCC_Combine_CTRP2,PCC_Combine_GDSC1,PCC_Combine_GDSC2)



Stats_Combine_CTRP1 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/Stats_Combine_CTRP1.csv",header=TRUE)
Stats_Combine_CTRP2 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/Stats_Combine_CTRP2.csv",header=TRUE)
Stats_Combine_GDSC1 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/Stats_Combine_GDSC1.csv",header=TRUE)
Stats_Combine_GDSC2 <- read.csv("C:/Users/EYC/Downloads/drive-download-20211127T042527Z-001/Stats_Combine_GDSC2.csv",header=TRUE)


Stats_Combine_CTRP1 <- melt(as.data.frame(Stats_Combine_CTRP1[-1]), id=c("method","drug"))
Stats_Combine_CTRP2 <- melt(as.data.frame(Stats_Combine_CTRP2[-1]), id=c("method","drug"))
Stats_Combine_GDSC1 <- melt(as.data.frame(Stats_Combine_GDSC1[-1]), id=c("method","drug"))
Stats_Combine_GDSC2 <- melt(as.data.frame(Stats_Combine_GDSC2[-1]), id=c("method","drug"))

Stats_Combine_CTRP1$dataset="CTRP1"
Stats_Combine_CTRP2$dataset="CTRP2"
Stats_Combine_GDSC1$dataset="GDSC1"
Stats_Combine_GDSC2$dataset="GDSC2"

stat <- rbind(Stats_Combine_CTRP1,Stats_Combine_CTRP2,Stats_Combine_GDSC1,Stats_Combine_GDSC2)
data = rbind(stat,pcc)
data <- data[data$value<=2 & data$value>=-0.3,] 


data <- data %>%
  mutate(flag = ifelse(method %in% c("rfG","rfinv1","rfinv2","rfinv3"), TRUE, FALSE),
         method_col = if_else(flag == TRUE, "1. RF family", "2. Others"))

vline <- data[!is.na(data$variable),] %>%
  group_by(variable) %>%
  summarise(mean=mean(value))


library(ggridges)
library(ggplot2)

data[!is.na(data$variable),] %>%
  arrange(value) %>%
  mutate(method = factor(method, levels=c("rfG","rfinv1","rfinv2",
                                          "rfinv3","rf","EN","ENglm","GR","KKNN","KNN","Lasso","pcr","pls","ranger","ridgeglm","svm","svmLinear2","treebag"))) %>%
ggplot(aes(x = value, y = method, fill = method,color=method_col)) +scale_y_discrete(limits=rev)+
  geom_density_ridges(quantile_lines = TRUE, alpha = 0.75,na.rm = TRUE,
                      quantiles = 2,size = 1)+facet_grid(~ variable,scales = "free")+#xlim(c(-1,3))
  geom_vline(data = vline, aes(xintercept = mean),linetype = "dashed",size=0.8)+
scale_color_manual(values = c("#F70020","#072C8F"))

## Summarize the stats by group 
(summ_bydataset <- data[!is.na(as.double(data$value)),] %>%
  group_by(dataset,variable,method) %>%
  summarise(stat.mean=mean(as.double(value)),stat.sd=sd(as.double(value))))

(summ_overall <- data[!is.na(as.double(data$value)),] %>%
    group_by(variable,method) %>%
    summarise(stat.mean=mean(as.double(value)),stat.sd=sd(as.double(value))))

reshape(summ_overall, idvar = "method", timevar = "variable", direction = "wide")

