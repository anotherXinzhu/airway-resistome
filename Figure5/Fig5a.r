

AMI <- data.table::fread("../_data/AMI.txt", data.table = F)
load("../_data/meta_80.RData")

plotD_full <- merge(meta, AMI, by.x = "SampleID", by.y = "Sample")



x = "A_FEV1_FVC_Post"  # A_FEV1_FVC_Post  in fig 5a，   O3_M_D10K in fig5c 
y = "total" # total in fig 5a, beta_lactam in fig 5c




library(dplyr)
plotD <- plotD_full %>%
  dplyr::select( SampleID, all_of(c(x,y)),  FEV1pred_post) %>%
  dplyr::rename(X = x, Y = y) 
head(plotD)
plotD <- plotD[complete.cases(plotD),]


# groups of population 
plotD$Group = sapply(1:nrow(plotD),
                     function(i){
                       if(plotD$X[i] < 70){
                         "COPD"
                       }else if(plotD$FEV1pred_post[i] < 80){
                         "PRISm"
                       }else{
                         "Healthy" 
                       }
                     })

library(ggplot2)
P_scatter <- ggplot( ) +
  geom_point(data=plotD, aes(x=X, y=Y, color=X), size=0.5) + 
  geom_smooth(data=plotD,aes(x=X, y=Y), method = "lm") +
  geom_smooth(data = plotD %>% filter(X > 70),
              aes(x=X, y=Y),
              color="red",
              method = "lm")+
  theme_bw() +
  theme(panel.grid = element_blank()) +
  xlab(x) + ylab(paste0(y)) +
  geom_vline(xintercept = 70, linetype="dashed") 


plotD_avg <- plotD %>%
  group_by(Group) %>%
  summarise(avg = mean(Y), sd=sd(Y), n=n())
P_avg.sd <- ggplot(data= plotD_avg, aes(x=Group, y=avg))  +
  geom_point(shape=23, size=3, color="orange", fill="orange") +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=0.1, color="orange") +
  theme_void() +
  theme(axis.text.y = element_text()) 




# stats :
wilcox.test(x = plotD$Y[plotD$Group == "COPD"],
            y = plotD$Y[plotD$Group != "COPD"])

wilcox.test(x = plotD$Y[plotD$Group == "PRISm"],
            y = plotD$Y[plotD$Group == "Healthy"])
