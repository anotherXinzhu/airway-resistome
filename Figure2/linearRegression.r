
library(data.table)
library(dplyr)


# ARG data  (take the total abundance and shannon as example) ---------
load("../_data/TotalARG.abund.shannon.RData")

totalARGIndices <- totalARGIndices %>%
  mutate(abund.log = log1p(ARGabundance)) %>%
  mutate(abund.st =( abund.log - mean(abund.log))/sd(abund.log))

head(totalARGIndices)
ARG.dat <- totalARGIndices %>% select(Sample, abund.st, ARGshannon) %>% rename(ARGAbundance = abund.st)

# meta data --------- 
load("../_data/meta_80.RData")
covars <- c("Age","BMI","Gender","District","medication")
covar.df <- meta.st %>% select(SampleID, all_of(covars))
meta.st <- meta.st %>% select(-all_of(covars))

# calculate linear regression between ARG features and exposure and health variables -----------------

Results <- NULL

for(arg.v in colnames(ARG.dat)[-1]){
  # arg.v = colnames(ARG.dat)[2]
  
  for(v in colnames(meta.st)[-1]){
    # v = "Smoking_binary"
    
    meta_c <- meta.st %>% select(SampleID,all_of(v))
    
    # 先去掉NA ----
    meta_c <- meta_c[complete.cases(meta_c),]
    if(length(unique(meta_c[,2])) == 1) next
    
    # 统一X和Y排序----
    ARG_c <- ARG.dat %>% select(Sample, all_of(arg.v))
    ARG_c <- ARG_c[match(meta_c$SampleID, ARG_c$Sample),]
    
    # statistics: lm without Covars ---------------
    testD <- cbind.data.frame(
      meta_c,
      ARGvalue = ARG_c[,2],
      stringsAsFactors=F
    ) %>% relocate(ARGvalue, .after = 'SampleID')
    
    
    # Lm with covars ---------------
    testD.covars <- merge(testD, covar.df, by = "SampleID")
    s1 <- lm(as.formula(paste0("ARGvalue ~ ", paste0(covars,collapse = "+" ), "+", colnames(testD)[3])), data = testD.covars)
    s0 <- lm(as.formula(paste0("ARGvalue ~ ", paste0(covars,collapse = "+" ))), data = testD.covars)
    an1 = anova(s1,s0)
    
    Lm.coef = s1$coef[colnames(testD)[3]]
    Lm.pval = an1[2,"Pr(>F)"]
    Lm.F = an1[2,"F"]
    Lm.r2 = summary(s1)$r.squared - summary(s0)$r.squared
    
    
    # bind results ------
    res_c <- data.frame(
      Y = arg.v,
      X = colnames(testD)[3], 
      Lm.coef,
      Lm.pval,
      Lm.F,
      Lm.r2
    )
    
    Results <- bind_rows(Results, res_c)
    
  }# loop through meta variables
}# loop through ARG features
