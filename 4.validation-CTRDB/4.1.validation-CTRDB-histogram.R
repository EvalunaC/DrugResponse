load("/extraspace/ychen42/validation-CTRDB/Epirubicin_testFrame_GDSC2.RData")
load("/extraspace/ychen42/validation-CTRDB/Epirubicin_trainFrame_GDSC2.RData")
load("/extraspace/ychen42/validation-CTRDB/GSE16446_median.RData")
load("/extraspace/ychen42/validation-CTRDB/GSE16446-processed.RData")#trainFrame_com, expr.mat,pt_pos,pt_neg

#write.csv(trainFrame,"/extraspace/ychen42/validation-CTRDB/Epirubicin_trainFrame_GDSC2.csv",quote = FALSE)
#write.csv(testFrame,"/extraspace/ychen42/validation-CTRDB/Epirubicin_testFrame_GDSC2.csv",quote = FALSE)


drug <- 'Epirubicin'


'''
run the models
get the predictions
'''

pred_DL <- read.csv("/extraspace/ychen42/validation-CTRDB/Epirubicin_testFrame_DLpredict.csv",header=F)

rownames(pred_DL) <- rownames(testFrame)
print(t.test(pred_DL[row.names(pred_DL)%in%rownames(pt_pos),], pred_DL[row.names(pred_DL)%in%rownames(pt_neg),]))

(wiltest <- wilcox.test(pred_DL[row.names(pred_DL)%in%rownames(pt_pos),], pred_DL[row.names(pred_DL)%in%rownames(pt_neg),]))
pdf("/extraspace/ychen42/validation-CTRDB/Epirubicin_DL_boxplot.pdf")
boxplot(list(Negative=pred_DL[row.names(pred_DL)%in%rownames(pt_pos),],  Positive=pred_DL[row.names(pred_DL)%in%rownames(pt_neg),]), las=1, col=c("#66c2a5", "#fc8d62"), pch=20, width=c(.75, .75), ylab="Predicted Epirubicin Sensitivity", xlab=paste("\n Wilcxn RS p = ", signif(wiltest$p.value,digits = 3), sep=""),cex.axis=.75, outcol="#00000033")
dev.off()

pt_pos <- rownames(pheno.anno[pheno.anno$`erbb2bimod:ch1`==1,])
pt_neg <- rownames(pheno.anno[pheno.anno$`erbb2bimod:ch1`==0,])

(wiltest <- wilcox.test(pred_DL[row.names(pred_DL)%in%pt_pos,], pred_DL[row.names(pred_DL)%in%pt_neg,]))
pdf("/extraspace/ychen42/validation-CTRDB/Epirubicin_DL_boxplot_erbb2bimod.pdf")
boxplot(list(Negative=pred_DL[row.names(pred_DL)%in%pt_pos,],  Positive=pred_DL[row.names(pred_DL)%in%pt_neg,]), las=1, col=c("#66c2a5", "#fc8d62"), pch=20, width=c(.75, .75), ylab="Predicted Epirubicin Sensitivity", xlab=paste("\n Wilcxn RS p = ", signif(wiltest$p.value,digits = 3), sep=""),cex.axis=.75, outcol="#00000033")
dev.off()



method_list <- c("GR paper linear Ridge",
#              "Random Forest",
              "Random Forest (Ranger)",
              "Principle Component Regression",
              "Partial Least Square",
              "Ridge GLM",
              "Lasso GLM",
#              "K-nearest neighbors ",
#              "Weighted K-nearest neighbors ",
              "Support Vector Machine",
#              "Treebag (bootstrap aggregating)",
              "Elastic Net",
              "RF+Lasso2019(out-of-bag)",
              "RF+Lasso2019(split-conformal)",
              "RF+Lasso2019(quantile regression forest)",
              "Neural Network"
              )

preds <- data.frame(preds_GR,
            preds_ranger$predictions,
            preds_pcr,
            preds_pls,
            preds_ridgeglm,
            preds_Lasso_1,
            preds_svm,
            preds_EN,
            model_rfinv1$testPred,
            model_rfinv2$testPred,
            model_rfinv3$testPred,
            pred_DL
          )
          colnames(preds) <- c("GR paper linear Ridge",
          #              "Random Forest",
                        "Random Forest (Ranger)",
                        "Principle Component Regression",
                        "Partial Least Square",
                        "Ridge GLM",
                        "Lasso GLM",
          #              "K-nearest neighbors ",
          #              "Weighted K-nearest neighbors ",
                        "Support Vector Machine",
          #              "Treebag (bootstrap aggregating)",
                        "Elastic Net",
                        "RF+Lasso2019(out-of-bag)",
                        "RF+Lasso2019(split-conformal)",
                        "RF+Lasso2019(quantile regression forest)",
                        "Neural Network"
                        )
preds$id <- rownames(preds)
preds_long <- melt(preds,id="id")
library(ggplot2)
preds_pos <- preds_long[preds_long$id%in%rownames(pt_pos),]
preds_neg <- preds_long[preds_long$id%in%rownames(pt_neg),]
preds_pos$HER2status <- "Positive"
preds_neg$HER2status <- "Negative"
preds_long <- rbind(preds_pos,preds_neg)

pdf("/extraspace/ychen42/validation-CTRDB/Epirubicin_split_violin.pdf", width=15, height=8)
ggplot(preds_long,
       aes(as.factor(variable), value, fill = HER2status),repel=TRUE) +
       scale_x_discrete(guide = guide_axis(angle = 20))+
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))
dev.off()


pdf("/extraspace/ychen42/validation-CTRDB/Epirubicin_boxplot.pdf", width=8, height=10)
par(mfrow=c(3,4))
for (i in 1:12){
pred <- as.data.frame(preds[[i]])
rownames(pred) <- rownames(testFrame)
(wiltest <- wilcox.test(pred[row.names(pred)%in%rownames(pt_pos),], pred[row.names(pred)%in%rownames(pt_neg),]))
vio(list(Negative=pred[row.names(pred)%in%rownames(pt_pos),],  Positive=pred[row.names(pred)%in%rownames(pt_neg),]), las=1, col=c("#66c2a5", "#fc8d62"), ylim=c(2.8,4.5), pch=20, width=c(.75, .75), ylab="Predicted Epirubicin Sensitivity", xlab=paste(method_list[i], "\n Wilcxn RS p = ", signif(wiltest$p.value,digits = 3), sep=""),cex.axis=.75, outcol="#00000033")
}
dev.off()

ggplot(diamonds[which(diamonds$cut %in% c("Fair", "Good")), ],
       aes(as.factor(color), carat, fill = cut)) +
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75))


(wiltest <- wilcox.test(pred[row.names(pred)%in%rownames(pt_pos),], sort(pred[row.names(pred)%in%rownames(pt_neg),])[1:57]))
