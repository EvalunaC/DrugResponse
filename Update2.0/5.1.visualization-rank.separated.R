facet_grid(~ Variety)


library(readxl)
library(ggplot2) # visualisation
library(reshape)


GDSC1.345drugs <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "GDSC1-345drugs")
GDSC2.192drugs <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "GDSC2-192drugs")
Total.GDSC <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "Total-GDSC")
CTRP1.354drugs <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "CTRP1-354drugs")
CTRP2.545drugs <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "CTRP2-545drugs")
CTRP <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "CTRP")
Total <- read_excel("C:/Users/EYC/iCloudDrive/Documents/GitHub/DrugResponse/stat_eva.xlsx",sheet  = "Total")

GDSC1.345drugs$dataset="GDSC1.345drugs"
GDSC2.192drugs$dataset="GDSC2.192drugs"
CTRP1.354drugs$dataset="CTRP1.354drugs"
CTRP2.545drugs$dataset="CTRP2.545drugs"

data <- rbind(GDSC1.345drugs,GDSC2.192drugs,CTRP1.354drugs,CTRP2.545drugs)


data_fmt <- data.frame(method=data$method,dataset=data$dataset,RMSE=-data$RMSE,
                       MAE=-data$MAE,R2=data$R2,PCC=data$`median of PCC`,PCC_0.5=data$`PCC>=0.5(#)`,AUC=data$AUC,AUC_median=data$`median of AUC`)



data_long <- melt(as.data.frame(data_fmt), id=c("method","dataset"))

library(dplyr)
data_long_rank <- data_long %>%
  group_by(variable,dataset) %>%
  mutate(rank = rank(-value))

data_long_rank2 <- data_long_rank %>%
  mutate(alpha = ifelse(method %in% c("rfG","rfinv1","rfinv2","rfinv3"), 0.9, 0.8),
         line_size = ifelse(method %in% c("rfG","rfinv1","rfinv2","rfinv3"), 2, 1.5))


ggplot(data = data_long_rank2, aes(x = variable, y = rank, group = method)) +
#  ggtitle("GDSC1 345drugs") +
  geom_line(aes(color = method, alpha = alpha), size = 2) +
  geom_point(aes(color = method, alpha = alpha), size = 4) +
  scale_y_reverse(breaks = 1:nrow(data_long_rank))+facet_grid(~ dataset)+scale_alpha(range = c(0.3, 9))


my_theme <- function() {
  # Colors
  color.background = "white"
  color.text = "#22211d"
  # Begin construction of chart
  theme_bw(base_size=15) +
    # Format background colors
    theme(panel.background = element_rect(fill=color.background, color=color.background)) +
    theme(plot.background  = element_rect(fill=color.background, color=color.background)) +
    theme(panel.border     = element_rect(color=color.background)) +
    theme(strip.background = element_rect(fill=color.background, color=color.background)) +
    # Format the grid
    theme(panel.grid.major.y = element_blank()) +
    theme(panel.grid.minor.y = element_blank()) +
    theme(axis.ticks       = element_blank()) +
    # Format the legend
    # Format title and axis labels
    theme(plot.title       = element_text(color=color.text, size=20, face = "bold")) +
    theme(axis.title.x     = element_text(size=14, color="black", face = "bold")) +
    theme(axis.title.y     = element_text(size=14, color="black", face = "bold", vjust=1.25)) +
    theme(axis.text.x      = element_text(size=10, vjust=0.5, hjust=0.5, color = color.text)) +
    theme(axis.text.y      = element_text(size=10, color = color.text)) +
    theme(strip.text       = element_text(face = "bold")) +
    # Plot margins
    theme(plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm"))
}


ggplot(data = data_long_rank2, aes(x = variable, y = rank, group = method)) +
  ggtitle("GDSC1 345drugs") +
  geom_line(aes(color = method_col, alpha = 0.9), size = 2) +
  geom_point(aes(color = method_col, alpha = 0.9), size = 4) +
  scale_y_reverse(breaks = 1:nrow(data_long_rank))+
  my_theme() +
  scale_color_manual(values = c("#F70020","#FB9701","#1A7D00","#072C8F","grey"))

#  scale_color_manual(values = c("#F70020","#191A1A","#FB9701","#1A7D00","#072C8F","grey"))


show.top.n <- 18
ggplot(data = data_long_rank2, aes(x = variable, y = rank, group = method)) +
  geom_line(aes(color = method_col, alpha = 1), size = 2) +
  geom_point(color = "#FFFFFF", size = 4) +
  geom_point(aes(color = method_col, alpha = 1), size = 4) +
  geom_point(color = "#FFFFFF", size = 1) +
  scale_y_reverse(breaks = 1:show.top.n) +
  scale_x_continuous(breaks = 1:18, minor_breaks = 1:18, expand = c(.05, .05)) +
  coord_cartesian(ylim = c(1,show.top.n)) + 
  theme(legend.position = "none") +
  labs(x = "Competition days with medals",
       y = "Rank",
       title = "PyeongChang 2018 Olympic Winter Games",
       subtitle = "Countries ranked by overall medals after each competition day") +
  my_theme() +
  scale_color_manual(values = c("#F70020","#191A1A","#FB9701","#1A7D00","#072C8F","grey"))

summary(data)

data[, as.list(summary(data)), by = method]
