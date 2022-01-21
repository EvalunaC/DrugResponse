source("Creating_drug_data.R")
library(gplots)
library("ggplot2")

load("/extraspace/ychen42/Drug_Response/VAEN/VAEN/Figure/Figure5/GSE32989/GSE32989.RData")
#======================== From 03.5B.R ===============================
EMT = read.table("/extraspace/ychen42/Drug_Response/VAEN/VAEN/Figure/Figure5/GSE32989/EMT.txt", as.is=T)
EMT = unique(EMT[,1])

match(EMT, rownames(expr.mat)) -> ii
ii = ii[!is.na(ii)]

prcomp(t(expr.mat[ii, ]))-> fit

which(fit$x[,1] < median(fit$x[,1])) -> idx
men.cells = rownames(fit$x)[idx]
epi.cells = rownames(fit$x)[-idx]
save(epi.cells, men.cells, file="/extraspace/ychen42/Drug_Response/VAEN/VAEN/Figure/Figure5/GSE32989/cells.RData")

pdf("/extraspace/ychen42/Drug_Response/yiqings_work/Method_validation_VAENdata/GSE32989_ZEB1_boxplot.pdf")
boxplot(expr.mat["ZEB1", epi.cells], expr.mat["ZEB1", men.cells])
dev.off()
######################
#pdf("/extraspace/ychen42/Drug_Response/yiqings_work/Method_validation_VAENdata/GSE32989_heatmap.pdf", #width=5, height=5)
#heatmap.2(expr.mat[ii, ], trace="none", col=greenred(100), cexRow=.5, cexCol=.5) -> h
#dev.off()
#cutree(as.hclust(h$colDendrogram),2) -> c1
########### ZEB1

dat = rbind( cbind( ZEB1=expr.mat["ZEB1", which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=expr.mat["ZEB1", which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

pdf("/extraspace/ychen42/Drug_Response/yiqings_work/Method_validation_VAENdata/GSE32989_ZEB1_boxplot.pdf")
p1 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() +
     labs(title=paste("GSE32989\n", "p = ", format(pvalue, digits=3)), x="", y = "ZEB1 expression") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1) )
print(p1)
dev.off()


#======================== Generate training and testing data ===============================
# Erlotinib, same with 1.3
#drug <- "Erlotinib"
#drug_data <- getDrugData(drug)
#trainFrame <- drug_data

val.RPKM = expr.mat

shared_genes = intersect(colnames(drug_data), rownames(expr.mat))
print( length(shared_genes) )
#[1] 11732

match(shared_genes, rownames(val.RPKM)) -> ii
new.RPKM = val.RPKM[ii, ]

#apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
#rownames(raw.RPKM) = shared_genes
#
#### scale to z, per-sample
#rank.GSE32989.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
#rank.GSE32989.gene.mat = apply(rank.GSE32989.gene.mat, 1, function(u){
#	u[which(u==1)] = 6162.5/6163
#	u
#})
#testFrame = apply(rank.GSE32989.gene.mat, 2, function(u){qnorm(u)} )

drug_data_nohead = drug_data[-1]
new.drug_data_nohead= drug_data_nohead[,match(shared_genes, colnames(drug_data_nohead))]
trainFrame <- data.frame(Resp = drug_data[,1], new.drug_data_nohead)
testFrame <- new.RPKM


#======================== HomogenizeData ===============================

trainingExprData=t(new.drug_data_nohead)
testExprData=new.RPKM
selection=1
batchCorrect="standardize"
removeLowVaryingGenes=0.2
removeLowVaringGenesFrom="rawData"
printOutput=TRUE


homData <- homogenizeData(testExprData, trainingExprData,
                          batchCorrect = batchCorrect, selection = selection, printOutput = printOutput)

trainFrame_homo <- data.frame(Resp = drug_data[,1], t(homData$train))
testFrame_homo <- data.frame(t(homData$test))





#======================== Get models from 1.1 ===============================





shared_genes = intersect(colnames(trainFrame), rownames(expr.mat))
print( length(shared_genes) )

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



homData <- homogenizeData(testExprData, trainingExprData,
                          batchCorrect = batchCorrect, selection = selection, printOutput = printOutput)
