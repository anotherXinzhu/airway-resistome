
AMI <- data.table::fread("../_data/AMI.txt", data.table = F)
load("../_data/meta_80.RData")

plotD_full <- merge(meta, AMI, by.x = "SampleID", by.y = "Sample")


x = "Biofuel_exposure"  # O3_M_D10K, Smoking_binary, Biofuel_exposure
y = "rifamycin" # beta_lactam, multidrug, rifamycin




library(dplyr)
plotD <- plotD_full %>%
  dplyr::select( SampleID, all_of(c(x,y)),  FEV1pred_post) %>%
  dplyr::rename(X = x, Y = y) 
head(plotD)
plotD <- plotD[complete.cases(plotD),]


library(ggplot2)

if(is.numeric(plotD$X)){
  P_scatter <- ggplot( ) +
    geom_point(data=plotD, aes(x=X, y=Y, color=X), size=0.5) + 
    geom_smooth(data=plotD,aes(x=X, y=Y), method = "lm") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab(x) + ylab(paste0(y)) +
    geom_vline(xintercept = median(plotD$X), linetype="dashed") 
  
  
}else{
  
  avg.dat <- 
    plotD %>% 
    group_by(X) %>% 
    summarise(avg = mean(Y), sd=sd(Y), n=n())  
  
  P_box <- ggplot(plotD, aes(x=X, y=Y)) +
    geom_boxplot(outliers = F, width =  0.6) +  #试了用X_2cat着色不好看
    geom_point(data=avg.dat, aes(x=X, y=avg),shape=23, size=3, color="orange", fill="orange") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    xlab(x) + ylab(paste0("mobile ",y," resistance genes %"))
  
  #stats：
  wilcox.test(Y ~ X, data = plotD)
}