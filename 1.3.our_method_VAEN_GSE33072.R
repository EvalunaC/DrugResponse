source("Creating_drug_data.R")

library(doMC)
options(cores = 35)
registerDoMC()
library("survminer")
library("survival")
#survival data
#From VAEN 5A.R
load("/extraspace/ychen42/Drug_Response/VAEN/VAEN/Figure/Figure5/GSE33072/GSE33072.RData")
grep("progression-free survival time",pheno.anno[,"characteristics_ch1.7"]) -> ii
pheno.anno[ii, ] -> pheno.anno.pfs
pts <- pheno.anno.pfs[,2]

#length(possibleDrugs2)

#match("Erlotinib",possibleDrugs2)
drug <- "Erlotinib"
drug_data <- getDrugData(drug)
trainFrame <- drug_data



#===================================================================================================
shared_genes = intersect(colnames(trainFrame), rownames(expr.mat))
print( length(shared_genes) )
#[1] 13163

new.expr.mat = t(expr.mat[match(shared_genes, rownames(expr.mat)),])
new.expr.mat2 = new.expr.mat[match(pts, rownames(new.expr.mat)),]

trainFrame_nohead = trainFrame[-1]
new.trainFrame_nohead= trainFrame_nohead[,match(shared_genes, colnames(trainFrame_nohead))]
new.trainFrame <- data.frame(Resp = trainFrame[,1], new.trainFrame_nohead)
trainFrame <- new.trainFrame
testFrame <- as.data.frame(new.expr.mat2)
#> dim(new.trainFrame)
#[1]   467 13164

#> trainFrame[1:5,1:5]
#          Resp        FTL     TMSB10       RPL19       RPS18
#A101D 12.81190 -0.1853226 -0.3161834 -0.25099430 -0.34407137
#A172  16.10201 -0.4273588  1.0818728 -0.57804263 -1.37570236
#A204  16.68315 -0.6823355 -0.4666072 -0.03047636  0.09378443
#A2058 12.88885  0.3754346 -0.8541696 -0.02088863  0.15786089
#A253  15.10721 -0.2188721  0.0731487  0.08990299  0.02970797

#> testFrame[1:5,1:5]
#               FTL   TMSB10    RPL19    RPS18    S100A6
#GSM677317 11.92248 13.48595 10.20855 12.80725 10.091832
#GSM677318 11.50687 13.39695 10.27374 12.64714  9.766183
#GSM677319 11.17122 13.46740 10.18908 12.78889  9.792214
#GSM677320 11.47420 13.04092 10.53985 12.63992  9.766870
#GSM677321 11.38385 13.42191 10.58061 12.54099 10.184466


#Prediction
preds_GR<-predict(model_GR,as.data.frame(testFrame))
preds_ranger<-predict(model_ranger,testFrame)
preds_pcr<-predict(model_pcr,testFrame)
preds_pls<-predict(model_pls,testFrame)
preds_ridgeglm <- predict(model_ridgeglm, s = best_lam_0, newx=as.matrix(testFrame))
preds_Lasso_1 <- predict(model_Lasso_1, s = best_lam_1, newx=as.matrix(testFrame))
preds_knn<-predict(model_KNN,as.data.frame(testFrame[1,]))
preds_svm<-predict(model_svm,testFrame)
preds_EN<-predict(model_EN,testFrame)
#preds_rf<-predict(model_rf,testFrame[1:2,])
preds_kknn<-predict(model_KKNN,testFrame)

preds<- list(preds_GR,
            preds_ranger,
            preds_pcr,
            preds_pls,
            preds_ridgeglm,
            preds_Lasso_1,
            preds_kknn,
            preds_svm,
            preds_EN)
methods <- c("GR paper linear Ridge",
            "Random Forest (Ranger)",
            "Principle Component Regression",
            "Partial Least Square",
            "Ridge GLM",
            "Lasso GLM",
            "Weighted KNN",
            "SVM",
            "Elastic Net")

#COX model
OS_YEAR = as.numeric(sapply(as.character(pheno.anno.pfs[,"characteristics_ch1.7"]), function(u)trimws(strsplit(u, split=":")[[1]][2])))
OS = sapply(as.character(pheno.anno.pfs[,"characteristics_ch1.6"]), function(u)trimws(strsplit(u, split=":")[[1]][2]))
Y1 = Surv(as.numeric(OS_YEAR), as.numeric(OS))

pdf(paste("/extraspace/ychen42/Drug_Response/yiqings_work/Method_validation_VAENdata/",drug,".pdf", sep=""))
cox <- list()
p <- list()
xvector <- list()
exp_coef <- list()
for (i in 1:9){
  xvector[[i]] = ifelse(preds[[i]] > mean(preds[[i]]), "HR", "LR")
  cox[[i]] <- summary(coxph(Y1 ~ preds[[i]]))
  p[[i]] <- cox[[i]]$coefficients[5]
  exp_coef[[i]] <- cox[[i]]$coefficients[2]
#  fit_binary = survfit(Surv(as.numeric(OS_YEAR), as.numeric(OS)) ~ xvector[[i]])
  dat = data.frame(cbind(OS_YEAR=OS_YEAR, OS=OS, X=xvector[[i]]))
  dat[,1] = as.numeric(as.character(dat[,1]))
  fit_binary = survfit( Surv(as.numeric(OS_YEAR), as.numeric(OS)) ~ xvector[[i]], data=dat)
  p_word <- paste(methods[[i]],"\n CoxReg w/ predicted response p =",format(p[[i]], digits=3),",\n exp(coef)=",format(exp_coef[[i]], digits=3), sep="")
  g1 <- ggsurvplot(fit_binary, data=dat,pval =p_word)     #, risk.table = TRUE,pval = TRUE,break.time.by = 1, ggtheme = theme_minimal())
  print(g1)
}
dev.off()
#===================================================================================================
