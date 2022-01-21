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


Total_fmt <- data.frame(method=Total$method,RMSE=-Total$RMSE,
                        MAE=-Total$MAE,R2=Total$R2,PCC=Total$`median of PCC`,AUC=Total$AUC,AUC_median=Total$`median of AUC`)

Total_long <- melt(as.data.frame(Total_fmt), id=c("method"))

library(dplyr)
Total_long_rank <- Total_long %>%
  group_by(variable) %>%
  mutate(rank = rank(-value))

ggplot(data = Total_long_rank, aes(x = variable, y = rank, group = method)) +
  ggtitle("Overall") +
  geom_line(aes(color = method, alpha = 1), size = 2) +
  geom_point(aes(color = method, alpha = 1), size = 4) +
  scale_y_reverse(breaks = 1:nrow(Total_long_rank))

Total_long_rank2 <- Total_long_rank %>%
  mutate(flag = ifelse(method %in% c("rfG","rfinv1","rfinv2","rfinv3"), TRUE, FALSE),
         method_col = if_else(flag == TRUE, method, "z_others"))

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

ggplot(data = Total_long_rank2, aes(x = variable, y = rank, group = method)) +
  ggtitle("Combining all four datasets") +
  geom_line(aes(color = method_col, alpha = 0.9), size = 2) +
  geom_point(aes(color = method_col, alpha = 0.9), size = 4) +
  scale_y_reverse(breaks = 1:nrow(Total_long_rank))+
  my_theme() +
  scale_color_manual(values = c("#F70020","#FB9701","#1A7D00","#072C8F","grey"))
#  scale_color_manual(values = c("#F70020","#191A1A","#FB9701","#1A7D00","#072C8F","grey"))

show.top.n <- 18
ggplot(data = Total_long_rank2, aes(x = variable, y = rank, group = method)) +
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
