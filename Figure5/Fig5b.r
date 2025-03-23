
AMI <- data.table::fread("../_data/AMI.txt", data.table = F)
load("../_data/meta_80.RData")

covars <- c("Age","BMI","Gender","District")
covar.df <- meta.st %>% select(SampleID, all_of(covars))

exposures <- c


# lm  ====================
Exposures <- c( "Smoking_binary","SHS_binary","Biofuel_exposure","SO2_M","NH4_M","CO_M_D10K", "Occupational_pollution","SO4_M_D10K", #这里CO和SOd4都留了
                'PM1_M',"year2pm25","year2pm10", # PM group都保留了
                "NO2_M_D10K", "NO3_M", "O3_M_D10K","CI_M") # NOx group只保留了NO2: O3
LungFunctions <- c("A_FEV1_FVC_Post","FEV1pred_post")
LungDisease <- c("Allergic_rhinitis","Rhinitis","Asthma","Tuberculosis","COPD_diagnosis","Chronic_Bronchitis")
LungSyptom <- c("CAT_score","Phlegm","Cough","Wheeze","Dyspnea")
OtherHealth <- c("Diabetes", "Osteoporosis", "CHD","Gastroesophageal_reflux","Anemia","Hypertension","Arrhythmia")
HealthCols <- c(LungFunctions, LungDisease, LungSyptom)


Results <- NULL


for(arg.v in colnames(AMI)[-1]){
  # arg.v = "total"
  
  for(v in c(Exposures, HealthCols)){
    # v = "Smoking_binary"
    
    meta_c <- meta.st %>% select(SampleID,all_of(v))
    
    # 先去掉NA ----
    meta_c <- meta_c[complete.cases(meta_c),]
    if(length(unique(meta_c[,2])) == 1) next
    
    # 统一X和Y排序----
    ARG_c <- AMI %>% select(Sample, all_of(arg.v))
    ARG_c <- ARG_c[match(meta_c$SampleID, ARG_c$Sample),]
    
    # statistics: lm  ---------------
    testD <- cbind.data.frame(
      meta_c,
      ARGvalue = ARG_c[,2],
      stringsAsFactors=F
    ) %>% relocate(ARGvalue, .after = 'SampleID')
    
    testD.covars <- merge(testD, covar.df, by = "SampleID")
    s1 <- lm(as.formula(paste0("ARGvalue ~ ", paste0(covars,collapse = "+" ), "+", colnames(testD)[3])), data = testD.covars)
    s0 <- lm(as.formula(paste0("ARGvalue ~ ", paste0(covars,collapse = "+" ))), data = testD.covars)
    an1 = anova(s1,s0)
    
    lm.coef = s1$coef[colnames(testD)[3]]
    lm.pval = an1[2,"Pr(>F)"]
    lm.r2 = summary(s1)$r.squared - summary(s0)$r.squared
    
    
    # bind results ------
    res_c <- data.frame(
      Y = arg.v,
      X = colnames(testD)[3], 
      lm.coef,
      lm.pval,
      lm.r2
    )
    
    Results <- bind_rows(Results, res_c)
    
  }# loop through meta variables
}# loop through ARG features


write.table(Results, file = "lm_AMI.meta.txt", sep = '\t', quote = F, row.names = F)


# plot the corr heatmaps ====================================
majorARGtypes <- 
  c("total","multidrug","peptide","tetracycline","quinolone",'aminocoumarin',"MLS","beta_lactam",
    "aminoglycoside","rifamycin","phosphonic acid", "diaminopyrimidine","phenicol") # with descending abundance from previous analysis 
  
  
# health variables ---------------
plotD.health <- Results %>%
  filter(X %in% c(HealthCols))  %>% 
  mutate(Zscore=-qnorm(lm.pval/2) *  sign(lm.coef)) %>% 
  filter(!is.na(Zscore)) %>%
  select(X, Y, lm.pval, Zscore, lm.coef) %>%
  mutate(fdr = p.adjust(lm.pval, method="fdr")) %>%
  mutate(sig.fdr = cut(fdr, 
                       breaks = c(0, 0.001, 0.01, 0.05, 1),
                       labels = c("***","**","*","")))

plotD.health$X <- factor(plotD.health$X, levels = HealthCols)
plotD.health$Y <- factor(plotD.health$Y, levels = rev(majorARGtypes))

ggplot(plotD.health) +
  geom_tile(aes(x=X, y=Y, fill=Zscore), color="#E5E5E5") +
  geom_text(aes(x=X, y=Y, label = sig.fdr), size=4, color="white") +
  scale_fill_gradientn(
    colors = c("blue", "#F7F7F7","red"),  # 自定义颜色
    values = scales::rescale(c(min(plotD.health$Zscore, na.rm = T), 
                               0, 
                               max(plotD.health$Zscore, na.rm = T))),  # 自定义值,
    #name = paste("CovarLm.coef\n","+/-: fdr < ",fdr.co),
    na.value = "white"   ) +   # 设置 NA 值的颜色
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))


# exposures ------------------
plotD.expos <-  Results %>%
  filter(X %in% Exposures)  %>% 
  mutate(Zscore=-qnorm(lm.pval/2) *  sign(lm.coef)) %>% 
  filter(!is.na(Zscore)) %>%
  select(X, Y, lm.pval, Zscore, lm.coef) %>%
  mutate(fdr = p.adjust(lm.pval, method="fdr")) %>%
  mutate(sig.fdr = cut(fdr, 
                       breaks = c(0, 0.001, 0.01, 0.05, 1),
                       labels = c("***","**","*","")))

dat.w <- plotD.expos %>% 
  reshape2::dcast(X~Y, value.var = "Zscore") %>% 
  tibble::column_to_rownames("X")

library(ggdendro)
if(T){
  df <- t(dat.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] 
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) 
  ddata_y <- dendro_data(dd.col) 
} 
order_env <- xx_names[[2]]
order_ARG <- xx_names[[1]]

plotD.expos$X <- factor(plotD.expos$X, levels = order_env)
plotD.expos$Y <- factor(plotD.expos$Y, levels = order_ARG)


library(ggpubr)

p0 <- ggplot(plotD.expos) +
  geom_tile(aes(x=X, y=Y, fill=Zscore), color="#E5E5E5") +
  geom_text(aes(x=X, y=Y, label = sig.fdr), size=3, color="white") +
  scale_fill_gradientn(
    colors = c("blue", "#F7F7F7","red"),  # 自定义颜色
    values = scales::rescale(c(min(plotD.expos$Zscore, na.rm = T), 
                               0, 
                               max(plotD.expos$Zscore, na.rm = T))),  # 自定义值,
    #name = paste("CovarLm.coef\n","+/-: fdr < ",fdr.co),
    na.value = "white"   ) +   # 设置 NA 值的颜色
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
  xlab("") + ylab("") 


p.expo <-  ggdendrogram(ddata_x, size = 2) + theme_dendro()
p1 <- ggarrange(ggplot(), p.expo, ggplot(), widths = c(0.15,0.55, 0.3),nrow = 1)

p.arg <-  ggdendrogram(ddata_y, rotate = TRUE, size = 2) + theme_dendro()
p.arg <- ggarrange(p.arg, ggplot(), nrow = 2, heights = c(0.7,0.3))
ggarrange(p1,
          ggarrange(p0, p.arg, widths = c(0.8, 0.2)),
          nrow = 2, heights = c(0.2, 0.8)) 
