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
