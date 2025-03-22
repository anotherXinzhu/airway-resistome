library(dplyr)

LungFunctions <- c("A_FEV1_FVC_Post","FEV1pred_post")
LungDisease <- c("Allergic_rhinitis","Rhinitis","Asthma","Tuberculosis","COPD_diagnosis","Chronic_Bronchitis")
LungSyptom <- c("CAT_score","Phlegm","Cough","Wheeze","Dyspnea")
HealthCols <- c(LungFunctions, LungDisease, LungSyptom)

pvals_meta.totalARGindex <- data.table::fread("Fig2a.sourceData_lm.txt", data.table = F)

plotD <- pvals_meta.totalARGindex %>% 
  filter(X %in% HealthCols)  %>%  
  mutate(zscore = -qnorm(lm.pval/2)) %>%
  mutate(sig =  cut(lm.pval,
                    breaks = c(-Inf,0.001, 0.01, 0.05,  Inf),
                    labels = c("***","**","*","")))


library(ggplot2)
plotD$X <- factor(plotD$X, levels = HealthCols)
plotD$zscore.dir <- plotD$zscore * sign(plotD$lm.coef)

ggplot(plotD) +
  geom_tile(aes(x=X, y=0, fill = zscore.dir), color="#E5E5E5") +
  scale_fill_gradientn(
    colors = c("#EB746A", "#F7F7F7","#47A9AE"),  # 自定义颜色
    values = scales::rescale(c(min(plotD$zscore.dir,na.rm = T), 
                               0, 
                               max(plotD$zscore.dir,na.rm = T))),  # 自定义值,
    name = "Correlation (zcore)",
    na.value = "white") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1) ) +
  ylab("") +
  facet_grid(Y~.)+
  geom_text(aes(x=X, y=0, label = sig))
