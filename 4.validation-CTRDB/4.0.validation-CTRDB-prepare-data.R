

setwd("/extraspace/ychen42/validation-CTRDB/")
load("/extraspace1/qli12/Data Source/GDSC2/GDSC2_Frame/Epirubicin_trainFrame.RData")
#load("C:/Users/EYC/Downloads/Epirubicin_trainFrame.RData")

library("GEOquery")
library("limma")

gse <- getGEO("GSE16446")

save(gse, file="GSE16446.RData")

exprs(gse[[1]]) -> datExpr0

########## section 1

library('preprocessCore')
normalize.quantiles(datExpr0) -> datExprLQ
dimnames(datExprLQ) = dimnames(datExpr0)

########## section 1 end


########## section 2

gpl = getGEO("GPL570")
Table(gpl) -> anno

match(rownames(datExprLQ), anno[,1]) -> ii
symbol = anno$`Gene Symbol`[ii]
apply(datExprLQ, 2, function(u)tapply(u, symbol, median)) -> expr.mat

com_gene <- intersect(colnames(trainFrame),rownames(expr.mat))

length(com_gene)/length(colnames(trainFrame))


pData(gse[[1]]) -> pheno.anno

save(expr.mat, pheno.anno, file="C:/Users/EYC/Downloads/GSE16446.RData")

table(pheno.anno$`her2fishbin:ch1`,pheno.anno$`erbb2bimod:ch1`)
pheno.anno$`erbb2bimod:ch1`

#     0  1
# 0  59  2
# 1   8 22
# NA 23  4

# HER2.bin: HER2 status by fluorescent in situ hybridization (FISH): 0:not amplified (ratio<2), 1: amplified (ratio>or=2)
pt_pos <- na.omit(pheno.anno[pheno.anno[,"her2fishbin:ch1"]==1 & pheno.anno[,"erbb2bimod:ch1"]==1,])
pt_neg <- na.omit(pheno.anno[pheno.anno[,"her2fishbin:ch1"]==0 & pheno.anno[,"erbb2bimod:ch1"]==0,])

trainFrame_com <- trainFrame[,com_gene]
expr.mat <- t(expr.mat[com_gene,])

#save(trainFrame_com, expr.mat,pt_pos,pt_neg, file="C:/Users/EYC/Downloads/GSE16446-processed.RData")


library(doMC)
options(cores = 45)
registerDoMC()


library(pRRophetic)
library(car)
library(ridge)

setwd("/extraspace/ychen42/validation-CTRDB/")
load("/extraspace/ychen42/validation-CTRDB/GSE16446-processed.RData")#trainFrame_com, expr.mat,pt_pos,pt_neg

homData <- homogenizeData(t(expr.mat), t(trainFrame_com),
                          batchCorrect = "standardize", selection = 1, printOutput = TRUE)

testFrame <- t(homData$test)
