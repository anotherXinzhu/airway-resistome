
vars2analyze <- c("Age", "BMI","Gender", "A_FEV1_FVC_Post","FEV1pred_post","CAT_score",
                  "Smoking_binary","Biofuel_exposure","Occupational_pollution","SHS_binary",
                  "year2pm25","year2pm10","CI_M","CO_M_D10K","NH4_M","NO2_M_D10K","NO3_M","O3_M_D10K",     
                  "PM1_M","SO2_M","SO4_M_D10K",   
                  "Allergic_rhinitis","Rhinitis","Asthma","Tuberculosis","Chronic_Bronchitis","COPD_diagnosis",
                  "Bronchiectasis", "IPF","Emphysema","OSAS",
                  "Phlegm","Cough","Wheeze","Dyspnea",
                  "Arrhythmia","Gastroesophageal_reflux","Hypertension","Anemia","Osteoporosis",
                  "HF","Depression","Diabetes","CHD","Stroke","Other_cancer","Medication")

load("../_data/ARGsubtype.abund.RData")
load("../_data/meta_80.RData")

library(dplyr)
ARGdat <- ARGsubtype.abund %>%
  tibble::column_to_rownames("ARGsubtype") %>%
  t() %>% as.data.frame()


library(vegan)

Adonis_res <- NULL
for(v in vars2analyze){
  
  dat.test <- meta %>% select(all_of(c("SampleID", v)))
  dat.test <- dat.test[complete.cases(dat.test),]
  
  common.sps <- intersect(rownames(ARGdat), dat.test$SampleID)
  
  ARGdat.tmp <- ARGdat[match(common.sps, rownames(ARGdat)),]
  dat.test <- dat.test[match(common.sps, dat.test$SampleID),]
  
  fml <- paste0("ARGdat.tmp ~ ", paste0(c(v), collapse = " + " ))
  Adonis <- try(adonis2(as.formula(fml), data = dat.test, permutations = 999, 
                        method="bray", by="margin"))
  
  if('try-error' %in% class(Adonis)){
    R2 <- NULL
    Fval <- NULL
    p <- NULL
  }else{
    R2 <- Adonis$R2[1]
    Fval <- Adonis$F[1]
    p <- Adonis$`Pr(>F)`[1]
  } 
  
  
  res <- c("variable" = v,
           "adonis.F" = Fval,
           "adonis.r2" = R2,
           "adonis.p" = p)
  
  Adonis_res <- bind_rows(Adonis_res, res)
}


write.table(Adonis_res, file = "AdonisResults.txt", sep = "\t", quote = F, row.names = F)


# plotting  -----------------------------
# 针对total ARG本身，以不同的variable类型着色：
vars.demogr <- c("Age","BMI","Gender" )

sigAdonis_res <-  Adonis_res %>% filter(fdr < 0.05) %>% filter(!variable %in% vars.demogr)
vars.exposure <- c("Smoking_binary","Biofuel_exposure","Occupational_pollution", "SHS_binary","year2pm25","year2pm10","CI_M","CO_M_D10K", 
                   "NH4_M","NO2_M_D10K","NO3_M","O3_M_D10K","PM1_M", "SO2_M","SO4_M_D10K")
vars.airwayhealth <- c("A_FEV1_FVC_Post","FEV1pred_post","CAT_score", "Allergic_rhinitis","Rhinitis","Asthma", "Tuberculosis","Chronic_Bronchitis",
                       "COPD_diagnosis","Bronchiectasis","IPF","Emphysema","OSAS","Phlegm","Cough","Wheeze","Dyspnea")
vars.lungFunc <-  c("A_FEV1_FVC_Post","FEV1pred_post","CAT_score")
vars.airwaySymp <- c("Phlegm","Cough","Wheeze","Dyspnea")
vars.airwayDiseas <- c("Allergic_rhinitis","Rhinitis","Asthma", "Tuberculosis","Chronic_Bronchitis", "COPD_diagnosis","Bronchiectasis","IPF","Emphysema","OSAS")
vars.otherDisease <- c( "Arrhythmia","Gastroesophageal_reflux", "Hypertension","Anemia","Osteoporosis","HF","Depression","Diabetes",
                        "CHD","Stroke","Other_cancer")

sigAdonis_res$var.type = 
  sapply(sigAdonis_res$variable, 
         function(x){
           if(x %in% vars.exposure){
             "Exposure"
           }else if(x %in% vars.lungFunc){
             "LungFunction"
           }else if(x %in% vars.airwayhealth){
             "AirwayHealthy"
           }else if(x %in% vars.otherDisease){
             "OtherDiseases"
           }
         })

sigAdonis_res$fdr.sig <-
  cut(sigAdonis_res$fdr,breaks = c(0,0.001,0.01,0.05,1), labels = c("***","**","*",""))

sigAdonis_res <- sigAdonis_res %>% arrange(adonis.F)
sigAdonis_res$variable <- factor(sigAdonis_res$variable, levels = sigAdonis_res$variable)

Colors <- setNames(c("#1D457FFF", "#61599DFF", "#C36377FF", "#EB7F54FF", "#F2AF4AFF"),
                   c("Demography","Exposure","LungFunction","AirwayHealthy", "OtherDiseases"))

ggplot(sigAdonis_res ) +
  geom_col(aes(x=adonis.F, y=variable, fill=var.type),
           position = "dodge2",  width = 0.8) +
  geom_text(aes(x=adonis.F, y=variable, label=fdr.sig),
            position = position_dodge(width = 0.9), hjust=-0.2) +
  scale_fill_manual(values = Colors) +
  theme_bw() + theme(panel.grid = element_blank())
