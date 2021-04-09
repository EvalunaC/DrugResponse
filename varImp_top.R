


top <- var_rf$importance[with(var_rf$importance,order(-Overall)),,drop = FALSE]
var10_rf <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),rf = seq(1,10))

top <- var_pcr$importance[with(var_pcr$importance,order(-Overall)),,drop = FALSE]
#var10_pcr <- top[1:10,,drop = FALSE]
var10_pcr <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),pcr = seq(1,10))

top <- var_pls$importance[with(var_pls$importance,order(-Overall)),,drop = FALSE]
var10_pls <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),pls = seq(1,10))

top <- var_ridgeglm[with(var_ridgeglm,order(-Overall)),,drop = FALSE]
var10_ridgeglm <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),ridgeglm = seq(1,10))

top <- var_Lasso_1[with(var_Lasso_1,order(-Overall)),,drop = FALSE]
var10_Lasso_1 <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),Lasso_1 = seq(1,10))

top <- var_KNN$importance[with(var_KNN$importance,order(-Overall)),,drop = FALSE]
var10_KNN <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),KNN = seq(1,10))

top <- var_svm$importance[with(var_svm$importance,order(-Overall)),,drop = FALSE]
var10_svm <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),svm = seq(1,10))

top <- var_treebag$importance[with(var_treebag$importance,order(-Overall)),,drop = FALSE]
var10_treebag <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),treebag = seq(1,10))

top <- var_EN$importance[with(var_EN$importance,order(-Overall)),,drop = FALSE]
var10_EN <- data.frame(gene=rownames(top[1:10,,drop = FALSE]),EN = seq(1,10))





data1 <- merge(var10_rf,var10_pcr,by="gene",all = T,suffixes = c("",""))
data2 <- merge(var10_pls,var10_ridgeglm,by="gene",all = T,suffixes = c("",""))
data3 <- merge(var10_Lasso_1,var10_KNN,by="gene",all = T,suffixes = c("",""))
data4 <- merge(var10_svm,var10_treebag,by="gene",all = T,suffixes = c("",""))

data5 <- merge(data1,data2,by="gene",suffixes = c("",""),all = T)
data6 <- merge(data3,data4,by="gene",suffixes = c("",""),all = T)
data7 <- merge(data5,data6,by="gene",suffixes = c("",""),all = T)
data <- merge(data7,var10_EN,by="gene",suffixes = c("",""),all = T)
#names(data)[names(data) == "Overall"] <- "Overall.EN"

data$overlap <- rowSums(!is.na(data))-1
data <- data[order(-data$overlap),]
