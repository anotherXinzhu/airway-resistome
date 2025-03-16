load("../_data/ARGsubtype.abund.RData")

ARG.avgAbund <- data.frame(ARGsubtype = ARGsubtype.abund$ARGsubtype,
                           avg.abund = rowMeans(ARGsubtype.abund[-1])) 

ARG.avgAbund <- merge(ARG.avgAbund, ARGinfo, by = 'ARGsubtype')

plotDat.subtypeNum = ARG.avgAbund %>%
  select(ARGsubtype, ResMechanism) %>%
  unique() %>%
  dplyr::group_by(ResMechanism) %>%
  dplyr::summarise(NumSubtype = n()) %>%
  as.data.frame() %>%
  arrange(desc(NumSubtype)) # 这里的Num指的是ARG subtype的数量
tp.lvls_byNum <- plotDat.subtypeNum$ResMechanism
plotDat.subtypeNum$ResMechanism <-  factor(plotDat.subtypeNum$ResMechanism, levels = tp.lvls_byNum)


plotDat.mechAbund = ARG.avgAbund %>%
  dplyr::group_by(ResMechanism) %>%
  dplyr::summarise(abund = sum(avg.abund)) %>%
  as.data.frame() %>%
  arrange(desc(abund))
plotDat.mechAbund$ResMechanism <- factor(plotDat.mechAbund$ResMechanism, levels = tp.lvls_byNum)

load("../_data/ARGresMech_colors.RData")


plotDat.mechAbund <- plotDat.mechAbund%>%
  arrange(ResMechanism) %>%
  mutate(fraction = abund/sum(abund)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  arrange(desc(abund)) 

P_abund <- ggplot(plotDat.mechAbund, 
                  aes(ymax=ymax, ymin=ymin, xmax=4, xmin=0, fill=ResMechanism)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  scale_fill_manual(values = resmech.colors) +
  coord_polar(theta="y") +
  xlim(c(0, 5)) + # 可调节donut thickness
  theme_void() +  
  theme(legend.position = "bottom") +
  ggtitle("Composition by ARG abundance")
P_abund


plotDat.subtypeNum <- plotDat.subtypeNum %>%
  arrange(ResMechanism) %>%
  mutate(fraction = NumSubtype/sum(NumSubtype)) %>%
  mutate(ymax= cumsum(fraction)) %>%
  mutate(ymin= c(0, head(ymax, n=-1))) %>%
  mutate(label_pos = (ymax + ymin) / 2) %>% 
  arrange(desc(NumSubtype))


P_num <- ggplot(plotDat.subtypeNum, 
                  aes(ymax=ymax, ymin=ymin, xmax=4, xmin=0, fill=ResMechanism)) + # xmin ~ xmax is where the ring is 
  geom_rect(color="white",lwd=0.5) +
  scale_fill_manual(values = resmech.colors) +
  coord_polar(theta="y") +
  xlim(c(0, 5)) + # 可调节donut thickness
  theme_void() +  
  theme(legend.position = "bottom") +
  ggtitle("Composition by ARG subtype number") 
P_num

