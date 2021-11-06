
query <- GDCquery(project = "TCGA-LGG",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")

GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

query <- GDCquery(project = "TCGA-ACC",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab",
                   file.type = "radiation")
GDCdownload(query)
clinical.BCRtab.radiation <- GDCprepare(query)

clinical.BCRtab.all$clinical_drug_acc  %>%
  head  %>%
 DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

 setwd("/extraspace/ychen42/Drug_Response/Data/")

 library(TCGAbiolinks)
 library(SummarizedExperiment)

 query <- GDCquery(project = "TCGA-LGG",
                   data.category = "Clinical",
                   data.type = "Clinical Supplement",
                   data.format = "BCR Biotab")
 GDCdownload(query)
 clinical.BCRtab.all <- GDCprepare(query)
 names(clinical.BCRtab.all)
 drug.lgg <- clinical.BCRtab.all$clinical_drug_lgg

 "death_days_to"
 clinical.BCRtab.all$clinical_drug_lgg$bcr_patient_uuid[1:5]


 clin.lgg <- GDCquery_clinic("TCGA-LGG", "clinical")

 pdf("/extraspace/ychen42/Drug_Response/Own_update2.0/test_TCGAanalyze_survival.pdf")
 TCGAanalyze_survival(clin.lgg,
                      "gender",
                      main = "TCGA Set\n LGG",height = 10, width=10)
 dev.off()

 lgg.outcome <- merge(clin,drug,by="id",all = FALSE)

 clin <- data.frame(id= clin.lgg$submitter_id,days_to_death = clin.lgg$days_to_death)
 drug <- data.frame(id= drug.lgg$bcr_patient_barcode,drug_name = drug.lgg$pharmaceutical_therapy_drug_name)
 data <-clin.lgg
 notDead <- is.na(data$days_to_death)
 data$s <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE)
 data <- data[, c("days_to_death", "s", "type")]
 if (any(notDead == TRUE)) {
        data[notDead, "days_to_death"] <- data[notDead, "days_to_last_follow_up"]
    }
    f.m <- formula(Surv(as.numeric(data$days_to_death), event = data$s) ~                                          data$type)
 fit <- do.call(survfit, list(formula = f.m, data = data))

 ##





 dataSurv <- TCGAanalyze_SurvivalKM(clinical_patient = dataClin,Survresult = TRUE)

 group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
 group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

 dataSurv <- TCGAanalyze_SurvivalKM(clinical_patient = dataClin,
                                    dataGE = dataFilt,
                                    Genelist = rownames(dataDEGs),
                                    Survresult = FALSE,
                                    ThreshTop = 0.67,
                                    ThreshDown = 0.33,
                                    p.cut = 0.05, group1, group2)

 clinical.BCRtab.radiation <- GDCprepare(query)

 clinical.BCRtab.all$clinical_drug_lgg  %>%
   head  %>%
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))






 library(SummarizedExperiment)
 library(TCGAbiolinks)
 query.exp <- GDCquery(project = "TCGA-LGG",
 legacy = TRUE,
 data.category = "Gene expression",
 data.type = "Gene expression quantification",
 platform = "Illumina HiSeq",
 file.type = "results",
 experimental.strategy = "RNA-Seq",
 #sample.type = c("Primary Tumor","Solid Tissue Normal"))
 sample.type = c("Primary solid Tumor"))
 GDCdownload(query.exp)
 lgg.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "lggExp.rda")

 # get subtype information
 dataSubt <- TCGAquery_subtype(tumor = "BRCA")

 # get clinical data
 dataClin <- GDCquery_clinic(project = "TCGA-BRCA","clinical")

 # Which samples are Primary Tumor
 dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP")
 # which samples are solid tissue normal
 dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")

 dataPrep <- TCGAanalyze_Preprocessing(object = brca.exp, cor.cut = 0.6)

 dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                       geneInfo = geneInfo,
                                       method = "gcContent")

 dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                   method = "quantile",
                                   qnt.cut =  0.25)

 dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                             mat2 = dataFilt[,dataSmTP],
                             Cond1type = "Normal",
                             Cond2type = "Tumor",
                             fdr.cut = 0.01 ,
                             logFC.cut = 1,
                             method = "glmLRT")
