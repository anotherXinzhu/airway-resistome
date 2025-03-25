ARGsubtype_abund <- data.table::fread("../_data/miceMetag_ARGabund.txt", data.table = F)
groups <- data.table::fread("../_data/miceExposure.txt",data.table = F)


# total abundance ---------------------------
ARGabund <- ARGsubtype_abund %>% summarise_if(is.numeric, sum, na.rm = TRUE) %>% t() %>% as.data.frame() %>% rename(ARGabund=V1)
plotD.abund <- merge(groups, ARGabund, by.x="Sample", by.y=0)
                  
plotD.abund$Group <- factor(plotD.abund$Group, levels = c( "Fresh Air","Biomass Fuel", "Cigarette Smoke"))

library(ggpubr)
ggboxplot(plotD.abund, x="Group", y="ARGabund", fill="Group", outlier.shape = NA) +
  geom_jitter(aes(x=Group, y=ARGabund,fill = Group), shape =21,  width = 0.25, alpha=0.6) +
  ylab("total ARG abundance (Cov/Gbp)") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_y_log10()+
  scale_fill_manual(values = c("gray","#E6949A","#47B1B6")) +
  scale_color_manual(values = c("gray","#E6949A","#47B1B6")) +
  stat_compare_means( aes(label = paste0("p = ", after_stat(p.format))),
                      ref.group = "Fresh Air",
                     method.args = list(alternative = "greater"),
                     label.y = log10(2000))


# AMI -------------------------------
AMI <- data.table::fread("../_data/miceMetag_AMI.txt", data.table = F)
plotD.AMI <- merge(groups, AMI, by="Sample")

plotD.AMI$Group <- factor(plotD.AMI$Group, levels = c( "Fresh Air","Biomass Fuel", "Cigarette Smoke"))
ggboxplot(plotD.AMI, x="Group", y="total", fill="Group", outlier.shape = NA) +
  geom_jitter(aes(x=Group, y=total,fill = Group), shape =21,  width = 0.25, alpha=0.6) +
  ylab("AMI of ARGs") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_y_log10()+
  scale_fill_manual(values = c("gray","#E6949A","#47B1B6")) +
  scale_color_manual(values = c("gray","#E6949A","#47B1B6")) +
  stat_compare_means( aes(label = paste0("p = ", after_stat(p.format))),
                      ref.group = "Fresh Air",
                      method.args = list(alternative = "greater"))


