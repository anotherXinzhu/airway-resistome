library(dplyr)

load("../_data/meta_80.RData")

covars <- c("Age","Gender","BMI")
covar_df <- meta.st %>% dplyr::select(SampleID, all_of(covars))

Exposures <- c( "NH4_M","CO_M_D10K","SO2_M", "Biofuel_exposure","Occupational_pollution","SHS_binary",
                "Smoking_binary","year2pm25",
                "O3_M_D10K","NO2_M_D10K") 
HealthCols <- c("A_FEV1_FVC_Post","FEV1pred_post", 
                "CAT_score","Phlegm","Cough" ,"Wheeze","Dyspnea", 
                "Allergic_rhinitis","Rhinitis","Asthma","Tuberculosis","Chronic_Bronchitis","COPD_diagnosis","Bronchiectasis","IPF","Emphysema","OSAS",
                "Arrhythmia","Gastroesophageal_reflux", "Hypertension","Anemia","Osteoporosis","HF","Depression","Diabetes", "CHD","Stroke","Other_cancer")
all(Exposures %in% colnames(meta.st))
all(HealthCols %in% colnames(meta.st))



# ARG data -----------------------------
load("../_data/ARGsubtype.abund.RData") # for the ARG information
load("../_data/ARGsubtype.abund_rmSingleton.RData")  # the ARGs with singleton removed
ARGDat <- ARGsubtype.abund 
ARGDat.st <- ARGsubtype.abund.st 
remove(ARGsubtype.abund, ARGsubtype.abund.st)

# modulation analysis =================
# Health ~ covars + Exposure * ARGsubtype
modulation_res <- NULL

for(healt in HealthCols){
  # healt = HealthCols[1]
  writeLines(paste0("Health var: ", which(HealthCols == healt),". ", healt, "============================"))
  
  for(exp in Exposures){
    # exp = Exposures[1]
    writeLines(paste0("     Exposure var: ", which(Exposures == exp),". ", exp, "--------------"))
    
    metaD <- meta.st[,c("SampleID", covars, exp, healt)]
    
    for(arg in colnames(ARGDat.st)){
      # arg = colnames(ARGDat.st)[1]
      
      argD <- ARGDat.st %>% dplyr::select(!!arg) %>% rename(ARG = !!arg)
      
      testD <- merge(metaD, argD, by.x = "SampleID", by.y = 0)
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



# significant ARG-Exposure interactions on Health indicators with ARG bringing worsening impact on Health indicators
modulation_res_sig <- 
  bind_rows(
  modulation_res %>% 
    filter(Health %in% c("FEV1pred_post","A_FEV1_FVC_Post")) %>%
    filter(Coefficient < 0 ),
  modulation_res %>%
    filter(!Health %in% c("FEV1pred_post","A_FEV1_FVC_Post")) %>%
    filter(Coefficient > 0 ))  %>% 
  filter(P.interaction < 0.05)  


# Additional criterias : ===========================================
# The population will be stratified by the median value of ARG abundance into ARG-high (ARG.H) and ARG-low (ARG.L) groups.
# If the exposure variable is categorical, Expo.H refers to the group with exposure, and Expo.L refers to the group without exposure.
# If the exposure variable is continuous, Expo.H refers to the high-exposure group, and Expo.L refers to the low-exposure group, stratified by the median value.

# 1) If both Exposure and Health are continuous variables:
#    The sample from ARG.H must show that the coefficient for exposure in lm(Health ~ Exposure) is > 0 and significant.
#    Additionally, the coefficient in ARG.H should be greater than that in ARG.L.
# 2) Otherwise: Exposure is treated as a categorical variable.
#   2.1) If Health is continuous, the wilcox.p in ARG.H must be significant,
#        and the health status in Expo.H should be lower than in Expo.L.
#   2.2) If Health is a categorical variable, the fisher.p in ARG.H must be significant,
#        and the Health% in Expo.H should be greater than in Expo.L.


SigRes <- modulation_res_sig

healthCat_criterias <- NULL
healthNum.expoCat_criterias <- NULL
healthNum.expoNum_criterias <- NULL
for(i in 1:nrow(SigRes)){
  # i=1
  
  arg = SigRes$ARG[i]
  expo = SigRes$Exposure[i]
  healt = SigRes$Health[i]
  
  # testD to determine the variable being categorical or continuous ------
  testD <- merge(meta %>% select(SampleID, !!expo, !!healt) %>% rename(Exposure=!!expo, Health=!!healt),
                 ARGDat %>% select(!!arg) %>% rename(ARG = !!arg),
                 by.x = "SampleID", by.y = 0)
  testD <- testD[complete.cases(testD),]
  testD$ARG.cat <- cut(testD$ARG,
                       breaks = c(-Inf, median(testD$ARG), Inf),
                       labels = c("ARG.L", "ARG.H"))
  testD.ARGL <- testD %>% filter(ARG.cat == 'ARG.L')
  testD.ARGH <- testD %>% filter(ARG.cat == "ARG.H")
  
  # testD.st to calculate interaction model --------
  testD.st <- merge(meta.st %>% select(SampleID, !!expo, !!healt) %>% rename(Exposure=!!expo, Health=!!healt),
                    ARGDat.st %>% select(!!arg) %>% rename(ARG = !!arg),
                    by.x = "SampleID", by.y = 0)
  testD.st <- testD.st[complete.cases(testD.st),]
  testD.st$ARG.cat <- cut(testD.st$ARG,
                          breaks = c(-Inf, median(testD.st$ARG), Inf),
                          labels = c("ARG.L", "ARG.H"))
  
  if(is.numeric(testD$Health)){ 
    testD.st.ARGL <- testD.st %>% filter(ARG.cat == 'ARG.L')
    testD.st.ARGH <- testD.st %>% filter(ARG.cat == "ARG.H")
    
    # Health is continous --------
    if(is.numeric(testD$Exposure)){
      # Exposure is continous
      lm.ARGL <- lm(Health ~ Exposure, data = testD.st.ARGL)
      lm.ARGH <- lm(Health ~ Exposure, data = testD.st.ARGH)
      ARGH.expoCoef <- summary(lm.ARGH)$coefficients['Exposure','Estimate']
      ARGH.expoPval <- summary(lm.ARGH)$coefficients['Exposure','Pr(>|t|)']
      ARGL.expoCoef <- summary(lm.ARGL)$coefficients['Exposure','Estimate']
      ARGL.expoPval <- summary(lm.ARGL)$coefficients['Exposure','Pr(>|t|)']
      
      healthNum.expoNum_res.c <- data.frame(
        ARGH.expoCoef,ARGH.expoPval,ARGL.expoCoef,ARGL.expoPval
      ) %>%
        mutate(Health = healt, Exposure=expo, ARG=arg,i.SigRes = i)
      healthNum.expoNum_criterias <- bind_rows(healthNum.expoNum_criterias, healthNum.expoNum_res.c)
    }else{
      if(length(unique(testD.st.ARGH$Exposure)) == 1) next
      
      m.wilcox = wilcox.test(Health ~ Exposure, data = testD.st.ARGH)
      if(healt %in% c("A_FEV1_FVC_Post", "FEV1pred_post","CAT_score")){
        healthNum.expoCat_res.c  <- testD.ARGH %>%
          group_by(Exposure) %>%
          summarise(avg.Health = mean(Health)) %>%
          reshape2::dcast(.~Exposure, value.var = 'avg.Health') %>%
          select(-`.`) %>%
          rename(ARGH.ExpoN.avgHealth = N,
                 ARGH.ExpoY.avgHealth = Y) %>%
          mutate(Health = healt, Exposure=expo, ARG=arg,i.SigRes = i) %>%
          mutate(ARGH.wilcoxPval = m.wilcox$p.value) 
        healthNum.expoCat_criterias <- bind_rows(healthNum.expoCat_criterias, healthNum.expoCat_res.c)
      }else{
        print("check Health variable：numeric but not lung function or CAT score")
        break
      }
      
    }
    
  }else{
    # Health is categorical --------
    
    if(is.numeric(testD$Exposure)){
      testD.st$Expo.cat <- cut(testD.st$Exposure,
                               breaks = c(-Inf, median(testD.st$Exposure), Inf),
                               labels = c("Expo.L", "Expo.H"))
    }else{
      testD.st$Expo.cat <- ifelse(testD.st$Exposure == 0, "Expo.L", "Expo.H")
    }
    
    testD.st.ARGL <- testD.st %>% filter(ARG.cat == 'ARG.L')
    testD.st.ARGH <- testD.st %>% filter(ARG.cat == "ARG.H")
    
    if(length(unique(testD.st.ARGH$Expo.cat)) == 1) next
    
    m.fisher <- fisher.test(table(testD.st.ARGH$Expo.cat, testD.st.ARGH$Health)) 
    m.fisher$p.value
    Health.perc <- testD.st.ARGH %>% 
      group_by(Expo.cat, Health) %>%
      summarise(n=n()) %>%
      mutate(total.n=sum(n)) %>%
      mutate(frac = n/total.n)  %>%
      reshape2::dcast(Health~Expo.cat, value.var = 'frac') %>%
      filter(Health == 1)
    Health.perc[is.na(Health.perc)] <- 0
    
    healthCat_res.c <- data.frame(
      ARGH.fisherPval = m.fisher$p.value,
      ARGH.ExpoH.HealtPerc = Health.perc$Expo.H,
      ARGH.ExpoL.HealtPerc = Health.perc$Expo.L
    ) %>% 
      mutate(Health = healt, Exposure=expo, ARG=arg, i.SigRes = i) 
    healthCat_criterias <- bind_rows(healthCat_criterias, healthCat_res.c)
  }
}



# Add a column to show whether to keep the link ---------------------------------------
# 1) If both Exposure and Health are continuous variables:
#    It is required that in the ARG.H sample, the coefficient for exposure in lm(Health ~ Exposure) must be > 0 and significant.
#    Additionally, the coefficient in ARG.H must be greater than the coefficient in ARG.L.
healthNum.expoNum_criterias$keep <-
  sapply(1:nrow(healthNum.expoNum_criterias),
         function(i){
           if(healthNum.expoNum_criterias$ARGH.expoPval[i] < 0.05 & 
              healthNum.expoNum_criterias$ARGH.expoCoef[i] > 0 &
              healthNum.expoNum_criterias$ARGH.expoCoef[i] > healthNum.expoNum_criterias$ARGL.expoCoef[i]){
             "keep"
           } else "discard"
         })

# If Exposure is categorical: 
# 2.1) If Health is continuous, it is required that the wilcox.p in ARG.H is significant,
#      and the health status in Expo.H should be lower than in Expo.L.healthNum.expoCat_criterias$keep <-
healthNum.expoCat_criterias$keep <-
  sapply(1:nrow(healthNum.expoCat_criterias),
         function(i){
           if(healthNum.expoCat_criterias$Health[i] %in% c("A_FEV1_FVC_Post","FEV1pred_post")){
             if(healthNum.expoCat_criterias$ARGH.wilcoxPval[i] < 0.05 &
                healthNum.expoCat_criterias$ARGH.ExpoY.avgHealth[i] < healthNum.expoCat_criterias$ARGH.ExpoN.avgHealth[i]){
               "keep"
             }else{
               "discard"
             }
           }else{
             if(healthNum.expoCat_criterias$ARGH.wilcoxPval[i] < 0.05 &
                healthNum.expoCat_criterias$ARGH.ExpoY.avgHealth[i] > healthNum.expoCat_criterias$ARGH.ExpoN.avgHealth[i]){
               "keep"
             }else{
               "discard"
             }
           }
         })

# 2.2) If Health is a categorical variable, it is required that the fisher.p in ARG.H is significant,
#      and the percentage of Health in Expo.H should be greater than in Expo.L.
healthCat_criterias$keep <- 
  sapply(1:nrow(healthCat_criterias),
         function(i){
           if(healthCat_criterias$ARGH.fisherPval[i] < 0.05 &
              healthCat_criterias$ARGH.ExpoH.HealtPerc[i] > healthCat_criterias$ARGH.ExpoL.HealtPerc[i]){
             "keep"
           }else{
             "discard"
           }
         })


Criterias <- bind_rows(healthCat_criterias %>% select(i.SigRes, Exposure, ARG, Health, keep),
                       healthNum.expoCat_criterias %>% select(i.SigRes, Exposure, ARG, Health, keep),
                       healthNum.expoNum_criterias %>% select(i.SigRes, Exposure, ARG, Health, keep))

SigRes <-
  merge(SigRes,
        Criterias %>% select(-i.SigRes),
        by=c("Exposure","ARG","Health")) %>% 
  filter(keep == "keep") %>% select(-keep)

OtherHealth <- c("Diabetes", "Osteoporosis", "CHD","Gastroesophageal_reflux","Anemia","Hypertension","Arrhythmia")
HealthCols <- HealthCols[!HealthCols %in% OtherHealth] # only take the airway health variables

SigRes <- SigRes %>% filter(Health %in% HealthCols) 
write.table(SigRes, 
            file = "significant_exposure.ARG_interactions.txt",
            sep = "\t", quote = F, row.names = F)

# plotting -------------------------
remove(list = ls())

LungFunctions <- c("A_FEV1_FVC_Post","FEV1pred_post")
LungDisease <- c("Allergic_rhinitis","Rhinitis","Asthma","Tuberculosis","COPD_diagnosis","Chronic_Bronchitis")
LungSyptom <- c("CAT_score","Phlegm","Cough","Wheeze","Dyspnea")
HealthCols <- c(LungFunctions, LungDisease, LungSyptom) # only take the airway health variables

SigRes <- data.table::fread("significant_exposure.ARG_interactions.txt",data.table = F)

load("../_data/ARGsubtype.abund.RData")
SigRes$ARGtype <- sapply(SigRes$ARG,
                         function(x){
                           ARGinfo$ARGtype[which(ARGinfo$ARGsubtype == x)]
                         })
SigRes$ARGtype %>% unique

load("../_data/ARGtype_colors.RData")

library(ggplot2)
library(ggpubr)

# plot by each exposure 
expo = "Smoking_binary"
plotD_hm <- SigRes %>% filter(Exposure %in% expo)
if(nrow(plotD_hm) == 0) next

ARGtypes <- plotD_hm %>% 
  select(ARG, ARGtype) %>%
  unique %>%
  arrange(ARG) %>%
  arrange(ARGtype)
ARGtypes$ARGtype[!ARGtypes$ARGtype %in% names(tp.colors)] <- "Others"

ARGimpact <- plotD_hm %>%
  group_by(ARG, ARGtype) %>%
  summarise(n.heal = n()) %>%
  arrange(desc(n.heal)) %>%
  arrange(ARGtype)

plotD_hm <- plotD_hm %>%
  reshape2::dcast(ARG ~ Health, value.var = "Coefficient") %>%
  reshape2::melt(id.vars = "ARG",variable.name = 'Health', value.name = "Coefficient")
plotD_hm$Health <- factor(plotD_hm$Health, levels = HealthCols)
plotD_hm$ARG <- factor(plotD_hm$ARG, levels = rev(ARGimpact$ARG))
ARGtypes$ARG <- factor(ARGtypes$ARG, levels = rev(ARGimpact$ARG))

p1 <- ggplot(plotD_hm) +
  geom_tile(aes(x=Health, y=ARG, fill=Coefficient), color="#E5E5E5") +
  scale_fill_gradientn(
    colors = c("blue", "#F7F7F7","red"),  
    values = scales::rescale(c(min(SigRes$Coefficient, na.rm = T), 
                               0, 
                               max(SigRes$Coefficient, na.rm = T))),  
    na.value = "white"  
  ) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5, color="black"),
        legend.title = element_text(size = 9),
        axis.text.y = element_text(size = 5, color="black")) 

p2 <- ggplot(ARGtypes) +
  geom_tile(aes(x=0, y=ARG, fill=ARGtype))+
  scale_fill_manual(values = tp.colors) +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())

ggarrange(p1, p2, widths = c(0.55, 0.45)) 

