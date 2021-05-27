library(ggplot2)
library("nparcomp")
library(fmsb)
library("dplyr")         

?mctp

#result0 <- read.csv("C:/Users/EYC/Downloads/Output_04302021_combine.csv", header = TRUE,na.strings="NA")
result0 <- read.csv("/Users/evaluna/Downloads/Output_04092021_combine2.csv", header = TRUE,na.strings="NA")
pcc <- read.csv("/Users/evaluna/Downloads/Output_cv28pcc_05132021_combine.csv", header = TRUE,na.strings="NA")

result <- merge(result0,pcc,by=c("method","drug"),all=TRUE)

result$method <- factor(result$method, levels = c("GR paper linear Ridge",
                                                  "Random Forest",
                                                  "Principle Component Regression",
                                                  "Partial Least Square",
                                                  "KNN",
                                                  "SVM",
                                                  "Elastic Net",
                                                  "Ridge GLM",
                                                  "Lasso GLM"))
p <- ggplot(result0, aes(factor(method), RMSE))

p <- p + geom_violin(aes(colour = "#1268FF"), alpha = 0.3)
q <- p + geom_violin(aes(y = R_Square, colour = "#3268FF"), alpha = 0.3)
q
old <- factor(result0$method, levels = c("virginica", "versicolor", "setosa"))
new <- levels(result0$method)
boxplot(RMSE~method,data=result0,ylim = c(0, 10),col=rainbow(10, alpha=0.2),angle = 75)

ggplot(result0, aes(x = RMSE)) +
  geom_histogram(fill = "white", colour = "black") +
  facet_grid(method ~ .)

anova(lm(RMSE ~ method , data = result0))
summary(lm(RMSE ~ method , data = result0))
plot(mctp(RMSE ~ method , data = result0))
plot(nparcomp(RMSE ~ method , data = result0))
summary(kruskal.test(RMSE ~ method , data = result0))
kruskal.test(RMSE ~ method , data = result0)$method

library("ggpubr")
result$RS <- result$R_Square



ggboxplot(result0, x = "method", y = "RMSE",
          color = "method", ylim = c(0, 10),
          ylab = "RMSE (Less is better)", xlab = "method",repel=TRUE)+  # outlier.shape = NA
  scale_x_discrete(guide = guide_axis(angle = 20))
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

ggboxplot(result, x = "method", y = "RS",
          color = "method", ylim = c(0, 1),
          ylab = "Simple R_square (Closer to 1 is better)", xlab = "method",repel=TRUE)+
  scale_x_discrete(guide = guide_axis(angle = 20))

ggboxplot(result0, x = "method", y = "R2",
          color = "method", ylim = c(0, 1),
          ylab = "R_square by correlation (Closer to 1 is better)", xlab = "method",repel=TRUE)+
  scale_x_discrete(guide = guide_axis(angle = 20))


ggboxplot(result0, x = "method", y = "MAE",
          color = "method", ylim = c(0, 5),
          ylab = "R_square (Less is better)", xlab = "method",repel=TRUE)+
  scale_x_discrete(guide = guide_axis(angle = 20))


ggboxplot(result, x = "method", y = "cv28_pcc",
          color = "method", ylim = c(0, 1),
          ylab = "PCC (Closer to 1 is better)", xlab = "method",repel=TRUE)+
  scale_x_discrete(guide = guide_axis(angle = 20))

tapply(result0$R_Square, result0$method, summary)



pcc <- read.csv("/Users/evaluna/Downloads/Output_cv28pcc_05132021_combine.csv", header = TRUE,na.strings="NA")
ggviolin(subset(pcc, !is.na(cv28_pcc)), x = "method", y = "cv28_pcc", 
          color = "method", ylim = c(0, 1),
          ylab = "PCC (Closer to 1 is better)", xlab = "method",
         repel=TRUE,draw_quantiles = 0.5,trim=TRUE,add = "jitter")+  # outlier.shape = NA
  scale_x_discrete(guide = guide_axis(angle = 20))+
  stat_summary(fun.data = function(x) data.frame(y=0.05, label =signif(mean(x),digits = 3)), geom="text") +
  theme(legend.position="none")
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

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
pcc=result_mean$pcc_mean
)


write.csv(result_mean,"/Users/evaluna/Downloads/result_mean.csv", quote=FALSE)
radarchart(result_per)


