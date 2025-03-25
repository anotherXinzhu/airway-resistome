
groups <- fread("../_data/miceExposure.txt", data.table = F)

# fold change of ARG subtypes
ARGsubtype.abund <- fread("../_data/miceMetag_ARGabund.txt",data.table = F)
CovGbp_ARG.l <- ARGsubtype.abund %>%
  reshape2::melt(id.vars = c( "ARGsubtype", "ARGtype"),
                 variable.name = "sample",
                 value.name = "ARGabund") %>%
  mutate(Group = sapply(sample,function(x) groups$Group[which(groups$Sample == x)])) 
head(CovGbp_ARG.l)

#calculate FC and pval for each ARGsubtype -------------

myComparisons <- 
  list(c("Biomass Fuel", "Fresh Air"),
       c("Cigarette Smoke", "Fresh Air"))

pvals_uniARGs <- NULL
for(cp in myComparisons){
  # cp = myComparisons[[1]]
  
  testD = CovGbp_ARG.l %>% 
    filter(Group %in% cp) 
  
  pval_c <- testD %>%
    group_by(ARGsubtype, ARGtype) %>%
    summarise(avg.Tm = mean(ARGabund[which(Group == cp[1])]),
              avg.Ct = mean(ARGabund[which(Group == cp[2])]),
              wilcox.p = wilcox.test(ARGabund ~ Group)$p.value,
              ttest.p = t.test(ARGabund ~ Group)$p.value) %>%
    ungroup() %>%
    mutate(FC = avg.Tm/avg.Ct) %>%
    mutate(Treatment=cp[1],
           Control=cp[2])
  
  pvals_uniARGs <- bind_rows(pvals_uniARGs, pval_c)
  
}


#  Reassign the control value of 0 to half of the minimum value, and recalculate the FC. 
ct.min = min(pvals_uniARGs$avg.Ct[pvals_uniARGs$avg.Ct != 0])
pvals_uniARGs$avg.Ct[pvals_uniARGs$avg.Ct == 0] <- ct.min/2
pvals_uniARGs$FC <- pvals_uniARGs$avg.Tm/pvals_uniARGs$avg.Ct


pvals_uniARGs <- pvals_uniARGs %>% 
  filter(FC != 0) %>%
  mutate(log2FC = log2(FC)) %>%
  filter(abs(log2FC) > log2(1.5)) %>%
  mutate(direction = ifelse(log2FC > 0, "increased","decreased")) 



# same pattern in the human respiratory dataset? 
expo.arg.cor_human <- fread("../Figure2/meta_ARGfamilies_lmRes.txt", data.table = F)
biofuelIncreased_insig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Biofuel_exposure_Y") %>%
     filter(lm.pval > 0.05) %>%
     filter(lm.coef > 0))$Y 

biofuelIncreased_sig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Biofuel_exposure_Y") %>%
     filter(lm.pval <= 0.05) %>%
     filter(lm.coef > 0))$Y 

biofuelDecreased_insig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Biofuel_exposure_Y") %>%
     filter(lm.pval > 0.05) %>%
     filter(lm.coef < 0))$Y 

biofuelDecreased_sig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Biofuel_exposure_Y") %>%
     filter(lm.pval <= 0.05) %>%
     filter(lm.coef < 0))$Y 

smokingIncreased_insig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Smoking_binary_Y") %>%
     filter(lm.pval > 0.05) %>%
     filter(lm.coef > 0) )$Y
smokingIncreased_sig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Smoking_binary_Y") %>%
     filter(lm.pval <= 0.05) %>%
     filter(lm.coef > 0) )$Y
smokingDecreased_insig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Smoking_binary_Y") %>%
     filter(lm.pval > 0.05) %>%
     filter(lm.coef < 0) )$Y
smokingDecreased_sig.human <- 
  (expo.arg.cor_human %>%
     filter(X == "Smoking_binary_Y") %>%
     filter(lm.pval <= 0.05) %>%
     filter(lm.coef < 0) )$Y



pvals_uniARGs$found.human = 
  sapply(1:nrow(pvals_uniARGs),
         function(i){
           if(pvals_uniARGs$Treatment[i] =="Biomass Fuel"){
             if(pvals_uniARGs$log2FC[i] > 0 & pvals_uniARGs$ARGsubtype[i] %in% biofuelIncreased_insig.human){
               "sameDirection inHuman"
             }else if(pvals_uniARGs$log2FC[i] < 0 & pvals_uniARGs$ARGsubtype[i] %in% biofuelDecreased_insig.human){
               "sameDirection inHuman"
             }else if(pvals_uniARGs$log2FC[i] > 0 & pvals_uniARGs$ARGsubtype[i] %in% biofuelIncreased_sig.human){
               "significant inHuman"
             }else if(pvals_uniARGs$log2FC[i] < 0 & pvals_uniARGs$ARGsubtype[i] %in% biofuelDecreased_sig.human){
               "significant inHuman"
             }else "not"
             
           }else if(pvals_uniARGs$Treatment[i] =="Cigarette Smoke"){
             if(pvals_uniARGs$log2FC[i] > 0 & pvals_uniARGs$ARGsubtype[i] %in% smokingIncreased_insig.human){
               "sameDirection inHuman"
             }else if(pvals_uniARGs$log2FC[i] < 0 & pvals_uniARGs$ARGsubtype[i] %in% smokingDecreased_insig.human){
               "sameDirection inHuman"
             }else if(pvals_uniARGs$log2FC[i] > 0 & pvals_uniARGs$ARGsubtype[i] %in% smokingIncreased_sig.human){
               "significant inHuman"
             }else if(pvals_uniARGs$log2FC[i] < 0 & pvals_uniARGs$ARGsubtype[i] %in% smokingDecreased_sig.human){
               "significant inHuman"
             }else "not"
           }
         })

# only plot the top 10 ARG types
load("ARGType_relAbundRanks.RData")
pvals_uniARGs <- pvals_uniARGs %>%
  filter(ARGtype %in% argTpRanks$ARGtype[1:10]) 

df_tile <- data.frame(x=unique(pvals_uniARGs$ARGtype), y=0)
df_bar <- pvals_uniARGs %>%
  group_by(ARGtype, Treatment)%>% 
  summarise(max = max(log2FC), 
            min=min(log2FC,0))


# plotting: overlapped with human, add label to the highest ones
ARG2show <- pvals_uniARGs %>%
  filter(direction=="increased") %>%
  filter(found.human %in% c("sameDirection inHuman","significant inHuman"))%>%
  group_by(ARGtype, Treatment) %>%
  slice_max(log2FC, n = 2)

load("../_data/ARGtype_colors.RData")

library(ggrepel)
ggplot() +
  geom_col(data = df_bar, aes(x = ARGtype,y = max),fill = "pink", alpha=0.5,width = 0.9) +
  geom_col(data = df_bar, aes(x = ARGtype,y = min),fill = "skyblue3", alpha=0.5,width = 0.9) +
  geom_jitter(data = pvals_uniARGs %>% filter(direction=="increased") %>% filter(found.human=="significant inHuman"), 
              aes(x = ARGtype, y = log2FC), 
              color="black", fill="red", stroke = 0.3,size = 1, 
              width =0.4, shape=21)+
  geom_jitter(data = pvals_uniARGs %>% filter(direction=="increased") %>% filter(found.human=="sameDirection inHuman"), 
              aes(x = ARGtype, y = log2FC), 
              color="red", fill="red",
              size = 0.7, width =0.4, shape=21)+
  geom_jitter(data = pvals_uniARGs %>% filter(direction=="increased") %>% filter(found.human=="not"), 
              aes(x = ARGtype, y = log2FC), 
              color="red", fill="white",
              size = 0.7, width =0.4, shape=21)+
  geom_jitter(data = pvals_uniARGs %>% filter(direction=="decreased") %>% filter(found.human=="significant inHuman"),
              aes(x = ARGtype, y = log2FC),
              color="black",fill="blue", stroke = 0.3, size = 1,
              width =0.4, shape=21)+
  geom_jitter(data = pvals_uniARGs %>% filter(direction=="decreased") %>% filter(found.human=="sameDirection inHuman"),
              aes(x = ARGtype, y = log2FC),
              color="blue",fill="blue",
              size = 0.7, width =0.4, shape=21)+
  geom_jitter(data = pvals_uniARGs %>% filter(direction=="decreased") %>% filter(found.human=="not"),
              aes(x = ARGtype, y = log2FC),
              color="blue",fill="white",
              size = 0.7, width =0.4, shape=21)+
  geom_tile(data = df_tile, aes(x = x,y = 0,fill = x), height = 0.5, show.legend = F) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1)) +
  #scale_color_manual(values = c("black","white")) +
  facet_grid(.~Treatment) +
  geom_text_repel(data = ARG2show, aes(x = ARGtype, y = log2FC, label = ARGsubtype),
                  box.padding = 0.35, max.overlaps = 10, size=2) +
  scale_fill_manual(values = tp.colors)
