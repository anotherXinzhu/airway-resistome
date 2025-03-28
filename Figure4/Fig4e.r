
load("../_data/HGTnetwork.RData")

# generate edges and nodes for network  ===============================

tmp1 <- HGT.edges %>%
  select(phylum, up.node) %>%
  mutate(down.node=up.node)%>%
  mutate(up.node=phylum)%>%
  select(down.node,up.node)%>%
  mutate(edgeType = "taxon")

plotD.edges <- bind_rows(
  tmp1,
  HGT.edges %>%
    select(up.node, down.node) %>%
    mutate(edgeType = "ARG")
) %>% unique


plotD.nodes <- bind_rows(
  cbind.data.frame(
    id = tmp1$up.node, 
    lvl = "phylum" 
  ) %>% unique,
  cbind.data.frame(
    id = tmp1$down.node,
    lvl = "tax.inHGT"
  ) %>% unique
)

test <- cbind.data.frame(
  id = HGT.edges$down.node,
  lvl = "ARG"
) %>% unique 

plotD.nodes <- bind_rows(plotD.nodes, test)


plotD.edges$up.node[!plotD.edges$up.node %in% plotD.nodes$id]
plotD.edges$down.node[!plotD.edges$down.node %in% plotD.nodes$id]


# coloring, sizing ----------------------
# edge color by phylum
plotD.edges$phylum <- 
  sapply(1:nrow(plotD.edges), 
         function(i){
           if(plotD.edges$edgeType[i] == "taxon"){
             plotD.edges$up.node[i]
           } else{
             HGT.edges$phylum[which(HGT.edges$up.node == plotD.edges$up.node[i])[1]]
           }
         })

# nodes color by genus
plotD.nodes <- merge(merge(plotD.nodes, HGT.edges %>% select(up.node, genus) %>% unique, by.x = "id", by.y = "up.node", all.x = T),
                     HGT.edges %>% select(up.node, phylum) %>% unique, 
                     by.x = "id", by.y = "up.node", all.x = T)


# in each phylum, assign color to a few top genus 
plotD.nodes$genus[which(plotD.nodes$lvl =="ARG")] <- "ARG"
plotD.nodes$phylum[which(plotD.nodes$lvl =="ARG")] <- "ARG"

test1 <- HGT.edges %>% group_by(phylum, genus) %>% summarise(n.connARG.genus=n()) # 每个phylum展示arg connection最多的几种genus，看看怎么设定阈值，比较好选颜色
test1 %>%
  group_by(phylum) %>%
  summarise(n.genus2show = sum(n.connARG.genus >= 5))%>%arrange(n.genus2show)
test2 <- HGT.edges %>% group_by(phylum, genus, species) %>% summarise(n.connARG.species=n())
test2 %>%
  group_by(phylum) %>%
  summarise(n.genus2show = sum(n.connARG.species >= 5))%>%arrange(n.genus2show)

plotD.nodes <- merge(merge(plotD.nodes, test1, by = c("phylum","genus"), all.x = T),
                     test2, by.x=c("id","phylum","genus"), by.y=c("species","phylum","genus"), all.x = T)


# For each phylum, only display the top 5 genera with the most ARG connections, and use the lightest color for the other minor ones.
test1 <- test1 %>%
  group_by(phylum) %>%
  mutate(colorBy = ifelse(n.connARG.genus > 5 & !grepl("unclassified", genus), genus, paste0("other.", unique(phylum)) )  ) %>%
  ungroup() %>%
  arrange(desc(n.connARG.genus)) %>% arrange(phylum)

test1$colorBy %>% unique()

load("../_data/HGTnetowrk.colors.RData")

plotD.nodes <- merge(plotD.nodes, test1, by = c("phylum", "genus","n.connARG.genus"), all.x = T)
plotD.nodes$colorBy[plotD.nodes$lvl == "ARG"] <- "ARG"
plotD.nodes$colorBy[plotD.nodes$lvl == "phylum"] <- plotD.nodes$id[plotD.nodes$lvl == "phylum"]
all(plotD.nodes$colorBy %in% names(all.colors))



# labels:
plotD.nodes$label<- NA
species2show <- test2$species[which(test2$n.connARG.species > 5)]
i.id <- which(plotD.nodes$lvl == "phylum" | plotD.nodes$id %in% species2show)
plotD.nodes$label[i.id] <- plotD.nodes$id[i.id]



# node size
tmp <- HGT.edges %>% 
  group_by(up.node) %>%
  summarise(n.conARG = n())
plotD.nodes <- merge(plotD.nodes, tmp %>% rename(nodeSize1 = n.conARG), by.x = "id", by.y="up.node", all.x = T)
plotD.nodes$nodeSize1[plotD.nodes$lvl == "phylum"] <-  max(tmp$n.conARG) * 1.2
plotD.nodes$nodeSize1[plotD.nodes$lvl == "ARG"] <-  min(tmp$n.conARG) * 0.5


# node shape: phylum: squared，tax.inHGT: circle; ARG: point 
plotD.nodes$lvl %>% unique

plotD.nodes <- plotD.nodes %>% relocate(id, .before = everything())


# build the network
library(igraph)
net <- graph_from_data_frame(d=plotD.edges, vertices=plotD.nodes, directed=F) 
net


library(ggraph)

ggraph(net, layout = 'igraph', algorithm="nicely") +
  geom_edge_link(aes(color=phylum,linetype=edgeType)) +
  geom_node_point( aes(size = nodeSize1, fill=colorBy, shape=lvl),color="black") + #可以设定size： numARG, log10(numARG)
  geom_node_text(aes(label = label), size=2, color="black", repel=T) +
  scale_edge_color_manual(values = color.phyla) +
  scale_fill_manual(values = all.colors) +
  # scale_color_manual(values = color.phyla) +
  scale_shape_manual(values = c(20,22,21)) +
  scale_size(range = c(2, 12)) + 
  theme_void()

# this is a rough plot for the HGT network, the refined layout are generated by Gephi and tailored on Adobe Illustrator.

