
#==================================================
#===   Evaluation stats boxplot
#==================================================

library("ggpubr")
library(ggplot2)
require(gridExtra)
library(cowplot)

Stat_table <- read.csv("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/192drugs_4Stats.csv", header = TRUE,na.strings="NA")

Stat_table$RMSE_negative <- -Stat_table$RMSE
Stat_table$MAE_negative <- -Stat_table$MAE

Stat_table<- Stat_table[!is.na(Stat_table$R2_corr),]

plot_RMSE_negative <- ggboxplot(Stat_table, x = "method", y = "RMSE_negative",
          color = "method", ylim = c(-5,0),font.label = list(size = 6),
          ylab = "Negative RMSE \n (Greater is better)")+ rremove("x.text")+ rremove("xlab")+ rremove("legend")
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))

plot_MAE_negative <- ggboxplot(Stat_table, x = "method", y = "MAE_negative",
          color = "method", ylim = c(-5,0),
          ylab = "Negative MAE \n(Greater is better)", xlab = "method",repel=TRUE)+ rremove("x.text")+ rremove("xlab")+ rremove("legend")

plot_R_Square <- ggboxplot(Stat_table, x = "method", y = "R2_simp",
          color = "method", ylim = c(0, 1),
          ylab = "Simple R2 \n(Closer to 1 is better)")+ rremove("x.text")+ rremove("xlab")+ rremove("legend")

plot_R2_corr <- ggboxplot(Stat_table, x = "method", y = "R2_corr",
          color = "method", ylim = c(0, 1),
          ylab = "R2 by correlation \n(Closer to 1 is better)", xlab = "method",repel=TRUE)+
  scale_x_discrete(guide = guide_axis(angle = 20))+ rremove("legend")

#plot_cv28_pcc <- ggboxplot(result, x = "method", y = "cv28_pcc",
#          color = "method", ylim = c(0.1, 0.8),font.label = list(size = 4),
#          ylab = "Pearson Cor Coef \n(Closer to 1 is better)", xlab = "method",repel=TRUE)+
#  scale_x_discrete(guide = guide_axis(angle = 20))+ rremove("legend")

#tapply(result0$R_Square, result0$method, summary)

filename = paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/192drugs_4Stats_boxplot07072021.pdf", sep="",collapse = NULL)
pdf(filename,width=10, height=15)
plot_grid(plot_RMSE_negative,
          plot_MAE_negative,
          plot_R_Square,
          plot_R2_corr,
          align = "v",
          nrow = 4,
          rel_heights = c(22/100,22/100,22/100,34/100))
dev.off()
#==================================================
#===   Potential:    1.0.4. Plotting prediction correlation (for one drug)
#==================================================
library("PerformanceAnalytics")
preds_df<- data.frame(GR=as.numeric(preds_GR),
            ridgeglm=as.numeric(preds_ridgeglm),
            Lasso=as.numeric(preds_Lasso_1),
            EN=as.numeric(preds_EN),
            rf=as.numeric(preds_rf),
            ranger=as.numeric(preds_ranger$predictions),
            rfinv1=as.numeric(model_rfinv1$testPred),
            rfinv2=as.numeric(model_rfinv2$testPred),
            rfinv3=as.numeric(model_rfinv3$testPred),
            pcr=as.numeric(preds_pcr),
            pls=as.numeric(preds_pls),
            knn=as.numeric(preds_knn),
            kknn=as.numeric(preds_kknn),
            svm=as.numeric(preds_svm),
            tree=as.numeric(preds_treebag)
          )
filename = paste("/extraspace/ychen42/Drug_Response/Own_update2.0/Evaluation/192drugs_preds_correlation.pdf", sep="",collapse = NULL)
pdf(filename,width=15, height=15)
chart.Correlation(preds_df, histogram=TRUE)
dev.off()
