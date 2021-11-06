
library(TCGAbiolinks)
library(SummarizedExperiment)

possibleDrugs2 <- read.csv("/extraspace/ychen42/Drug_Response/Data/192Drug_list.csv")


project <- c("TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
"TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM")

surv_result <- data.frame()

#TCGA-LGG:
for (l in 3:length(project)){
#cate <- "TCGA-LGG"
cate <- project[l]
clin <- GDCquery_clinic(cate, "clinical")

query <- GDCquery(project = cate,
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab",
                legacy=TRUE)

GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")
names(clinical.drug)
drug <- data.frame(submitter_id= clinical.drug$bcr_patient_barcode,drug_name = clinical.drug$drug_name)

clinic_data <- merge(clin,drug,by="submitter_id")
dim(clinic_data)
a <- table(clinic_data$drug_name)
drug_names <- intersect(names(a[a>200]),as.vector(possibleDrugs2$x))#More than 50 pts.
print(a[a>200])
print(drug_names)
#TCGA-ACC     [1] 18 73     pass
#TCGA-BLCA    [1] 303  76
