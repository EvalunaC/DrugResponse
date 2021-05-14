library(ggplot2)
library("nparcomp")
?mctp

result0 <- read.csv("C:/Users/EYC/Downloads/Output_04302021_combine.csv", header = TRUE,na.strings="NA")
result0 <- read.csv("/Users/evaluna/Output_04092021_combine.csv", header = TRUE,na.strings="NA")
p <- ggplot(result0, aes(factor(method), RMSE))
p <- p + geom_violin(aes(colour = "#1268FF"), alpha = 0.3)
q <- p + geom_violin(aes(y = R_Square, colour = "#3268FF"), alpha = 0.3)
q


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
result0$RS <- result0$R_Square



ggboxplot(result0, x = "method", y = "RMSE",
          color = "method", ylim = c(0, 10),
          ylab = "RMSE (Less is better)", xlab = "method",repel=TRUE)+  # outlier.shape = NA
  scale_x_discrete(guide = guide_axis(angle = 20))
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

ggboxplot(result0, x = "method", y = "RS",
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

tapply(result0$R_Square, result0$method, summary)
