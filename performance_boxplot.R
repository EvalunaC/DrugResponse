library(ggplot2)
library("nparcomp")
library(fmsb)
library("dplyr")         

?mctp

#result0 <- read.csv("C:/Users/EYC/Downloads/Output_cv28_05272021.csv", header = TRUE,na.strings="NA")
result0 <- read.csv("/Users/evaluna/Downloads/cv28_05272021.csv", header = TRUE,na.strings="NA")
pcc <- read.csv("/Users/evaluna/Downloads/Output_cv28pcc_05132021_combine2.csv", header = TRUE,na.strings="NA")



result <- merge(result0,pcc,by=c("methods","drug"),all=TRUE)
result <- result[!is.na(result$RMSE),]
result <- result[!is.na(result$cv28_pcc),]
#result <- result[result$drug!="AZD8055",]

result$methods <- factor(result$methods, levels = c("GR paper linear Ridge",
                                                  "Random Forest",
                                                  "Principle Component Regression",
                                                  "Partial Least Square",
                                                  "KNN",
                                                  "SVM",
                                                  "Elastic Net",
                                                  "Ridge GLM",
                                                  "Lasso GLM"))

result <- subset(result,select=-c(X.1.x,X.x,X.1.y,X.y))
summary(result)

result$RMSE_negative <- -result$RMSE
result$MAE_negative <- -result$MAE



ggplot(result, aes(x = RMSE)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(method ~ .)

anova(lm(RMSE ~ method , data = result0))
summary(lm(RMSE ~ method , data = result0))
plot(mctp(RMSE ~ method , data = result0))
plot(nparcomp(RMSE ~ method , data = result0))
summary(kruskal.test(RMSE ~ method , data = result0))
kruskal.test(RMSE ~ method , data = result0)$method

library("ggpubr")


plot_RMSE_negative <- ggboxplot(result, x = "methods", y = "RMSE_negative",
          color = "methods", ylim = c(-5,0),font.label = list(size = 6),
          ylab = "Negative RMSE \n (Greater is better)")+ rremove("x.text")+ rremove("xlab")+ rremove("legend")
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

plot_R_Square <- ggboxplot(result, x = "methods", y = "R_Square",
          color = "methods", ylim = c(-0.3, 0.45),
          ylab = "Simple R2 \n(Closer to 1 is better)")+ rremove("x.text")+ rremove("xlab")+ rremove("legend")

plot_R2_corr <- ggboxplot(result, x = "methods", y = "R2",
          color = "methods", ylim = c(0, 0.6),
          ylab = "R2 by correlation \n(Closer to 1 is better)")+ rremove("x.text")+ rremove("xlab")+ rremove("legend")


plot_MAE_negative <- ggboxplot(result, x = "methods", y = "MAE_negative",
          color = "methods", ylim = c(-4,0),
          ylab = "Negative MAE \n(Greater is better)", xlab = "methods",repel=TRUE)+ rremove("x.text")+ rremove("xlab")+ rremove("legend")


plot_cv28_pcc <- ggboxplot(result, x = "methods", y = "cv28_pcc",
          color = "methods", ylim = c(0.1, 0.8),font.label = list(size = 4),
          ylab = "Pearson Cor Coef \n(Closer to 1 is better)", xlab = "method",repel=TRUE)+
  scale_x_discrete(guide = guide_axis(angle = 20))+ rremove("legend")

#tapply(result0$R_Square, result0$method, summary)

require(gridExtra)
library(ggplot2)
library(cowplot)
pdf("/Users/evaluna/Downloads/evaluation_result/Output_cv28_05272021.pdf",width=7,height=10)
#grid.arrange(plot_RMSE_negative,
#             plot_R_Square,
#             plot_R2_corr,
#             plot_MAE_negative,
#             plot_cv28_pcc, 
#             ncol=1,nrow=5)

plot_grid(plot_RMSE_negative,
          plot_R_Square,
          plot_R2_corr,
          plot_MAE_negative,
          plot_cv28_pcc, 
          align = "v", 
          nrow = 5, 
          rel_heights = c(1/6, 1/6,1/6,1/6,1/3))
dev.off()


pcc <- read.csv("/Users/evaluna/Downloads/Output_cv28pcc_05132021_combine.csv", header = TRUE,na.strings="NA")
ggviolin(subset(pcc, !is.na(cv28_pcc)), x = "methods", y = "cv28_pcc", 
          color = "methods", ylim = c(0, 1),
          ylab = "PCC (Closer to 1 is better)", xlab = "method",
         repel=TRUE,draw_quantiles = 0.5,trim=TRUE,add = "jitter")+  # outlier.shape = NA
  scale_x_discrete(guide = guide_axis(angle = 20))+
  stat_summary(fun.data = function(x) data.frame(y=0.05, label =signif(median(x),digits = 3)), geom="text") +
  theme(legend.position="none")
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))
ggboxplot(subset(pcc, !is.na(cv28_pcc)), x = "methods", y = "cv28_pcc", 
         color = "methods", ylim = c(0, 1),
         ylab = "PCC (Closer to 1 is better)", xlab = "methods with median values labeled",
         repel=TRUE,draw_quantiles = 0.5,trim=TRUE,add = "jitter")+  # outlier.shape = NA
  scale_x_discrete(guide = guide_axis(angle = 20))+
  stat_summary(fun.data = function(x) data.frame(y=0.05, label =signif(median(x),digits = 3)), geom="text") +
  theme(legend.position="none")



radarchart(result)

create_beautiful_radarchart(
  data = result, caxislabels = c(0, 5, 10, 15, 20),
  color = c("#00AFBB", "#E7B800", "#FC4E07")
)

result_clean <- data.frame()

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}


result_mean <- result %>%
  group_by(method) %>%
  summarize(RMSE_mean = mean(RMSE),mae_mean = mean(MAE),r2simple_mean = mean(RS))

r2cor_mean <- subset(result, !is.na(R2)) %>%
  group_by(method) %>%
  summarize(r2cor_mean = mean(R2))

pcc_mean <- subset(result, !is.na(cv28_pcc)) %>%
  group_by(method) %>%
  summarize(pcc_mean = mean(cv28_pcc))

result_mean$r2cor_mean <- r2cor_mean$r2cor_mean
result_mean$pcc_mean <- pcc_mean$pcc_mean

result_per <- data.frame(method=result_mean$method,rmse=(max(result_mean$RMSE_mean)-result_mean$RMSE_mean/(-min(result_mean$RMSE_mea),
mae=(1-(result_mean$mae_mean/(max(result_mean$mae_mean)-min(result_mean$mae_mean)))),
r2simple =result_mean$r2simple_mean,
r2cor=result_mean$r2cor_mean,
pcc=result_mean$pcc_mean)


write.csv(result_mean,"/Users/evaluna/Downloads/result_mean.csv", quote=FALSE)
radarchart(result_per)

