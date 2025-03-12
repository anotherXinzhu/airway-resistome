load("../_data/ARGsubtype.abund.RData")


ARG.avgAbund <- data.frame(ARGsubtype = ARGsubtype.abund$ARGsubtype,
                           avg.abund = rowMeans(ARGsubtype.abund[-1])) 

ARG.avgAbund <- merge(ARG.avgAbund, ARGinfo, by = 'ARGsubtype')

plotDat.subtypeNum = ARG.avgAbund %>%
  select(ARGtype, ARGsubtype) %>%
  unique() %>%
  dplyr::group_by(ARGtype) %>%
  dplyr::summarise(NumSubtype = n()) %>%
  as.data.frame() %>%
  arrange(desc(NumSubtype)) # 这里的Num指的是ARG subtype的数量
top10subtype <- plotDat.subtypeNum$ARGtype[1:10]



plotDat.typeAbund = ARG.avgAbund %>%
  dplyr::group_by(ARGtype) %>%
  dplyr::summarise(abund = sum(avg.abund)) %>%
  as.data.frame() %>%
  arrange(desc(abund))
top10Abund <- plotDat.typeAbund$ARGtype[1:10]

load("../_data/ARGtype_colors.RData")


plotDat.typeAbund$ARGtype[!plotDat.typeAbund$ARGtype %in% c(top10Abund,top10subtype) ] <- "Others"
plotDat.typeAbund <- plotDat.typeAbund %>%
  dplyr::group_by(ARGtype) %>%
  dplyr::summarise(abund = sum(abund)) %>%
  as.data.frame() %>% 
  arrange(desc(abund)) 
tp.lvls <- c(plotDat.typeAbund$ARGtype[plotDat.typeAbund$ARGtype != "Others"],"Others")


plotDat.subtypeNum$ARGtype[!plotDat.subtypeNum$ARGtype %in% c(top10Abund,top10subtype) ] <- "Others"
plotDat.subtypeNum <- plotDat.subtypeNum %>%
  dplyr::group_by(ARGtype) %>%
  dplyr::summarise(NumSubtype = sum(NumSubtype)) %>%
  as.data.frame() %>% 
  arrange(desc(NumSubtype)) 
tp.lvls_byNum <- c(plotDat.subtypeNum$ARGtype[plotDat.subtypeNum$ARGtype != "Others"],"Others")


plotDat.typeAbund$ARGtype <- factor(plotDat.typeAbund$ARGtype, levels = tp.lvls_byNum)  

plotDat.typeAbund <- plotDat.typeAbund%>%
  arrange(ARGtype) %>%
  mutate(fraction = abund/sum(abund)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  arrange(desc(abund)) 

P_abund <- ggplot(plotDat.typeAbund, 
                    aes(ymax=ymax, ymin=ymin, xmax=4, xmin=0, fill=ARGtype)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  scale_fill_manual(values = tp.colors) +
  coord_polar(theta="y") +
  xlim(c(0, 5)) + # 可调节donut thickness
  theme_void() +  
  theme(legend.position = "bottom") +
  ggtitle("ARG composition by ARG abundance")




plotDat.subtypeNum$ARGtype <- factor(plotDat.subtypeNum$ARGtype, levels = tp.lvls_byNum)   

plotDat.subtypeNum <- plotDat.subtypeNum %>%
  arrange(ARGtype) %>%
  mutate(fraction = NumSubtype/sum(NumSubtype)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  arrange(desc(NumSubtype))

P_num <- ggplot(plotDat.subtypeNum, 
                          aes(ymax=ymax, ymin=ymin, xmax=4, xmin=0, fill=ARGtype)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  scale_fill_manual(values = tp.colors) +
  coord_polar(theta="y") +
  xlim(c(0, 5)) + # 可调节donut thickness
  theme_void() +  
  theme(legend.position = "bottom") +
  ggtitle("ARG composition by ARG subtype number") 


