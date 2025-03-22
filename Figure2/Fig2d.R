



pvals_meta.totalARGindex <- data.table::fread("Fig2a.sourceData_lm.txt", data.table = F)
load("exposure_correlations_plotdat.RData")

Exposures <- c( "Smoking_binary","SHS_binary","Biofuel_exposure","SO2_M","NH4_M","CO_M_D10K", "Occupational_pollution","SO4_M_D10K", 
                'PM1_M',"year2pm25","year2pm10",
                "NO2_M_D10K", "NO3_M", "O3_M_D10K","CI_M") 


Colors = c("#47B1B6","#E6949A")
names(Colors) <- c("1", "-1")  

assocLimits <- c(0.04, 5.2)  #  minimum and maximum value of association
lineWid.range <- c(0.1, 1.3)  # line width corresponding to the minimum and maximum value of association 
alpha.range <- c(0.2, 1)  # alpha value corresponding to the minimum and maximum value of association 


plotD_curves <- cbind.data.frame(
  link1.start = "ARG abundance",  # link1 represents the association between ARG abundance and all exposures
  link1.end = levels(corr_exposures$Var1),
  link1.x1 = 3/5*length(levels(corr_exposures$Var1)), # X coordinate of ARG abundance
  link1.x2 = seq(1,length(levels(corr_exposures$Var1)),1), 
  link1.y1 = length(levels(corr_exposures$Var1)) - 1,  #  Y coordinate of ARG abundance
  link1.y2 = seq(length(levels(corr_exposures$Var1)),1,-1), 
  
  link2.start = "ARG shannon",  # link1 represents the association between ARG shannon and all exposures
  link2.end = levels(corr_exposures$Var1),
  link2.x1 =  length(levels(corr_exposures$Var1)) - 1,  #  X coordinate of ARG shannon
  link2.x2 = seq(1,length(levels(corr_exposures$Var1)),1), 
  link2.y1 =  3/5*length(levels(corr_exposures$Var1)),  # Y coordinate of ARG shannon
  link2.y2 = seq(length(levels(corr_exposures$Var1)),1,-1)
  
) 




link1.info <- pvals_meta.totalARGindex %>% 
  filter(Y=="ARGAbundance") %>%
  filter(X %in% Exposures) %>%
  select(X, lm.coef, lm.pval) %>%
  rename(link1.coeff = lm.coef, 
         link1.pval = lm.pval) %>%
  mutate(link1.zscore =  -qnorm(link1.pval/2) ) %>%
  mutate(link1.fdr = p.adjust(link1.pval, method="fdr")) %>%
  mutate(link1.fdr.sig = ifelse(link1.fdr < 0.05, "*", "ns"))

link2.info <- pvals_meta.totalARGindex %>% 
  filter(Y=="ARGshannon") %>% 
  select(X, lm.coef, lm.pval) %>%
  rename(link2.coeff = lm.coef, 
         link2.pval = lm.pval) %>%
  mutate(link2.zscore =  -qnorm(link2.pval/2) )%>%
  mutate(link2.fdr = p.adjust(link2.pval, method="fdr")) %>%
  mutate(link2.fdr.sig = ifelse(link2.fdr < 0.05, "*", "ns"))



plotD_curves <- 
  merge(merge(plotD_curves,
              link1.info,
              by.x="link1.end", by.y="X"), 
        link2.info, 
        by.x = "link2.end", by.y = "X")

plotD_curves$link1.end <- factor(plotD_curves$link1.end, levels = levels(corr_exposures$Var1))
plotD_curves <- plotD_curves %>% arrange(link1.end)



link1.curve.max = 0.2
link1.curve.min = 0.1
link2.curve.max= -0.1
link2.curve.min = -0.2

plotD_curves$link1.curvation <- seq(link1.curve.max,link1.curve.min,length.out=length(levels(corr_exposures$Var1)))
plotD_curves$link2.curvation <- seq(link2.curve.max,link2.curve.min,length.out=length(levels(corr_exposures$Var1)))


tmp1 <- plotD_curves %>% select(contains("link1")) %>% mutate(type = "link1")
tmp2 <- plotD_curves %>% select(contains("link2")) %>% mutate(type = "link2")
colnames(tmp1) <- sub("link1.","",colnames(tmp1))
colnames(tmp2) <- sub("link2.","",colnames(tmp2))
plotD_curves.l <- bind_rows(tmp1, tmp2)



# Plot the curves one by one and assign a different curvature value to each curve.
library(ggplot2)
curves_hostt <- lapply(1:nrow(plotD_curves.l), function(i) {
  geom_curve(aes(x = plotD_curves.l$x1[i], y = plotD_curves.l$y1[i],
                 xend = plotD_curves.l$x2[i], yend = plotD_curves.l$y2[i],  
                 color = as.character(sign(plotD_curves.l$coeff[i])),
                 linewidth = plotD_curves.l$zscore[i], 
                 alpha=plotD_curves.l$zscore[i]),  
             curvature = plotD_curves.l$curvation[i]) 
})



# plot the correlation heatmap between exposure factors
p1_hm <- ggplot(corr_exposures) +
  geom_tile(aes(x=var1.lvl, y=var2.lvl, fill =lm.Coeff ), color="white") +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_fill_gradientn(
    colors = c("#542788", "#E5E5E5","#B35806"),  
    values = scales::rescale(c(min(corr_exposures$lm.Coeff,na.rm = T), 
                               0, 
                               max(corr_exposures$lm.Coeff,na.rm = T))),  # 自定义值,
    #name = paste("CovarLm.coef\n","+/-: fdr < ",fdr.co),
    na.value = "white"   ) +
  geom_text(aes(x=var1.lvl, y=var2.lvl, label = fdr.txt))


# add correlation curves to the heatmap 
p_curves <- p1_hm + 
  curves_hostt +
  scale_x_continuous(breaks = 0:16,  
                     labels = c("", levels(corr_exposures$Var1), "") # manually add x axis text
  ) + 
  scale_y_continuous(breaks = 0:16,  
                     labels = c("", rev(levels(corr_exposures$Var1)), "")  # manually add y axis text
  ) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text.y = element_text()) +
  scale_linewidth(limits = assocLimits, range = lineWid.range) +  
  scale_alpha(limits = assocLimits, range = alpha.range) + 
  scale_color_manual(values = Colors) +
  guides(alpha = guide_legend(title = "Correlation (zscore)"),
         linewidth = guide_legend(title = "Correlation (zscore)"),
         color = guide_legend(title = "Correlation"))


plotD_text <- plotD_curves.l %>%
  select(start, x1, y1) %>%
  unique
p_curves <- p_curves +
  geom_text(data = plotD_text, aes(x=x1,y=y1, label = start)) 
p_curves

