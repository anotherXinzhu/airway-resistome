library(dplyr)
load("../_data/meta_80.RData")

LungFunctions <- c("A_FEV1_FVC_Post","FEV1pred_post")
LungDisease <- c("Allergic_rhinitis","Rhinitis","Asthma","Tuberculosis","COPD_diagnosis","Chronic_Bronchitis")
LungSyptom <- c("CAT_score","Phlegm","Cough","Wheeze","Dyspnea")
OtherHealth <- c("Diabetes", "Osteoporosis", "CHD","Gastroesophageal_reflux","Anemia","Hypertension","Arrhythmia")

HealthCols <- c(LungFunctions, LungDisease, LungSyptom)


Exposures <- c("NH4_M","CO_M_D10K","SO2_M", "Biofuel_exposure","Occupational_pollution","SHS_binary","SO4_M_D10K", 
               "Smoking_binary",'PM1_M',"year2pm25","year2pm10",
               "NO2_M_D10K", "NO3_M", "O3_M_D10K","CI_M") 


# ARG data -------------------------------------
ARG.dat <- data.table::fread("../_data/AMI.txt",data.table=F)

# modulation analysis =================
# Health ~ covars + Exposure * transPerc
covars <- c("Age","Gender","BMI")
covar_df <- meta.st %>% dplyr::select(SampleID, all_of(covars))

modulation_res <- NULL

for(healt in HealthCols){
  # healt = HealthCols[1]
  writeLines(paste0("Health var: ", which(HealthCols == healt),". ", healt, "============================"))
  
  for(exp in Exposures){
    # exp=Exposures[1]
    writeLines(paste0("     Exposure var: ", which(Exposures == exp),". ", exp, "--------------"))
    
    metaD <- meta.st[,c("SampleID", covars, exp, healt)]
    
    for(arg in colnames(ARG.dat)[-1]){
      # arg = colnames(ARG.dat)[-1][1]
      argD <- ARG.dat %>% dplyr::select(Sample, !!arg) %>% rename(ARG = !!arg)
      
      testD <- merge(metaD, argD, by.x = "SampleID", by.y = "Sample")
      testD <- testD[complete.cases(testD),] 
      
      fml <- paste0(healt, " ~ ", paste(covars, collapse = " + "), " + ",
                    exp, " + ", "ARG", " + ",
                    exp, " * ", "ARG" );fml
      m <- lm(as.formula(fml), data = testD)
      i.coef = grep(":", rownames( summary(m)$coefficients), fixed = T)
      if(length(i.coef) == 0 ) next
      coeff = summary(m)$coefficients[i.coef,"Estimate"]
      
      an1 = anova(m)
      i.an1 = grep(":",rownames(an1), fixed = T)
      F.stat = an1[i.an1,"F value"] 
      Pvalue = an1[i.an1,"Pr(>F)"]  
      
      res_c <- data.frame("Health" = healt,
                          "Exposure" = exp,
                          "ARG" = arg, 
                          "Coefficient" = coeff,
                          "F.value" = F.stat,
                          "P.interaction" = Pvalue)
      modulation_res <- bind_rows(modulation_res, res_c)
    }# loop through ARGs
  }# loop through Exposures
}# loop through Health cols


#  ARG-Exposure interactions on Health indicators with ARG bringing worsening impact on Health indicators
modulation_res$ARG.effect <- 
  sapply(1:nrow(modulation_res),
         function(i){
           if(modulation_res$Health[i] %in% c("FEV1pred_post","A_FEV1_FVC_Post")){
             if(modulation_res$Coefficient[i] < 0){
               "worseningHealth"
             }else {
               "relieving"
             }
           }else{
             if(modulation_res$Coefficient[i] > 0){
               "worseningHealth"
             }else {
               "relieving"
             }
           }
         })


SigRes <- modulation_res %>% 
  filter(P.interaction < 0.1) %>%
  filter(ARG.effect == "worseningHealth") %>%
  select(-ARG.effct)

write.table(SigRes, 
            file = "significant_exposure.AMI_interactions.txt",
            sep = "\t", quote = F, row.names = F)

# plotting ==============================
remove(list = ls())
SigRes <- data.table::fread("significant_exposure.AMI_interactions.txt",data.table = F)


SigRes$Health <- factor(SigRes$Health, 
                        levels =c("A_FEV1_FVC_Post", "FEV1pred_post", "Allergic_rhinitis", "COPD_diagnosis","Phlegm", "Cough", "Wheeze", "Dyspnea"))

# only AMI of two ARG types have significant interaction effect with exposure: 
tp = "total"  # total or MLS

plotD_sub <- SigRes %>% filter(ARG==tp)
ggplot(plotD_sub) +
  geom_tile(aes(x=Health, y=Exposure, fill=Coefficient), color="#E5E5E5", linewidth=0.35) +
  # geom_text(aes(x=Health, y=Y, label=label.fdr), color="black", size=3)+
  scale_fill_gradientn(
    colors = c("blue", "#F7F7F7","red"),  # 自定义颜色
    values = scales::rescale(c(min(SigRes$Coefficient, na.rm = T), 
                               0, 
                               max(SigRes$Coefficient, na.rm = T))),  # 自定义值,
    #name = paste("CovarLm.coef\n","+/-: fdr < ",fdr.co),
    na.value = "white"  # 设置 NA 值的颜色
  ) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 5, color="black"),
        legend.title = element_text(size = 9),
        axis.text.y = element_text(size = 5, color="black"),
        legend.position = "none") 
