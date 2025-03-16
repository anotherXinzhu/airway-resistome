
library(dplyr)

load("../_data/ARG.host.combined.RData")
taxlvls <- colnames(ARG.host.combined.update)[2:8] 

for(tl in taxlvls){
  # tl = taxlvls[1]
  
  ARGabund.df <- ARG.host.combined.update %>%
    # select(all_of(tl), uniACC, avgGbp) %>%
    group_by(!!sym(tl)) %>%  # group_by的column名存在tl中，可以不用先去重新定义column名，直接这样使用
    summarise(ARGavgGbp = sum(avgGbp))
  
  assign(paste0("numARG.",tl), ARGabund.df, envir = .GlobalEnv)
}

# Find the corresponding relationships for different taxonomic levels
taxLinks <- ARG.host.combined.update %>%
  select(all_of(taxlvls)) %>%
  unique

# edges ------------------------------------------
edges.dat <- NULL
for(ti in 1:6){
  # ti = 1
  ti.n = ti +1
  tl = taxlvls[ti]
  tl.n = taxlvls[ti.n]
  
  edges <- taxLinks[,c(tl, tl.n)] %>% unique
  colnames(edges) <- c("up.node","down.node")
  
  up.numARG <- eval(parse(text = paste0("numARG.",tl)))
  colnames(up.numARG) <- c("up.node", "up.ARGavgGbp")
  
  down.numARG <- eval(parse(text = paste0("numARG.",tl.n)))
  colnames(down.numARG) <- c("down.node",  "down.ARGavgGbp")
  
  edges.dat_c <- 
    merge(merge(edges, down.numARG, by="down.node", all=T),
          up.numARG, by="up.node", all=T) %>%
    mutate(up.lvl = tl, down.lvl=tl.n) 
  
  edges.dat <- bind_rows(edges.dat, edges.dat_c)
}


edges.dat$up.lvl <- factor(edges.dat$up.lvl, levels = taxlvls[1:6])
edges.dat<- edges.dat %>% arrange(desc(down.ARGavgGbp)) %>% arrange(up.node) %>% arrange(up.lvl)

taxlvl.cols <- c( "#BA7999FF","#EFC86EFF", "#97C684FF", "#6F9969FF", "#AAB5D5FF", "#5C66A8FF", "#E69B99FF", "black")
names(taxlvl.cols) <- c(taxlvls,"ARG")

taxlvl.size <- c(7,6,5,4,3,2,1, 1)
names(taxlvl.size) <- c(taxlvls,"ARG")


# For each upper taxonomy level, remove the lower taxonomy levels that contribute less than 1%
plotD.edges <- edges.dat%>%
  filter(!is.na(up.node)) %>% filter(!is.na(down.node)) %>%
  dplyr::group_by(up.node) %>%
  dplyr::mutate(frac =  down.ARGavgGbp / sum(down.ARGavgGbp) ) %>% 
  filter(frac >= 0.01) %>%  
  ungroup() %>%
  arrange(desc(down.ARGavgGbp)) %>% arrange(up.node) %>% arrange(up.lvl)



# nodes data---------------------------------------------------
library(scales)
plotD.nodes <- rbind(
  plotD.edges %>% 
    mutate(id=up.node, lvl=up.lvl, numARG = up.ARGavgGbp) %>%
    select(id, lvl, numARG) %>%
    as.data.frame(),
  plotD.edges %>% 
    mutate(id=down.node, lvl=down.lvl, numARG = down.ARGavgGbp) %>%
    select(id, lvl, numARG) %>%
    as.data.frame()
) %>% unique %>%
  mutate(lvl.col = sapply(lvl, function(x) taxlvl.cols[x])) %>% # for plot()
  mutate(lvl.size = sapply(lvl, function(x) taxlvl.size[x])) %>% # for plot()
  mutate(numARG.size = rescale(log10(numARG), to = c(1, 8)))  # for plot()
head(plotD.nodes)



# label topN genus with most ARGs
# Among the top fifteen genera, the 10th genus "Bulleidia" and the 12th genus "Pauljensenia" belong to phyla with contributions less than 1%, so they will not be displayed in the network diagram.
# thus extend to 16th and 17th genera
topNGenus <- 17
topGenus <- (plotD.nodes %>%
               filter(lvl=="genus") %>% filter(!grepl("g_unclassified", id)) %>%
               arrange(desc(numARG)))$id[1:topNGenus]



# color unclassified taxonomy  with gray
plotD.nodes <- plotD.nodes %>% filter(id %in% c(plotD.edges$up.node, plotD.edges$down.node))
plotD.nodes$lvl<-factor(plotD.nodes$lvl,levels = c("domain","phylum","class","order","family","genus","species","Unclassified"))
plotD.nodes$lvl[grep("_unclassified", plotD.nodes$id)] <- "Unclassified"

taxlvl.cols <- c(taxlvl.cols, "gray")
names(taxlvl.cols)[length(taxlvl.cols)] <- "Unclassified"


# node labels
plotD.nodes$labels <- sapply(1:nrow(plotD.nodes),
                            function(i){
                              if(plotD.nodes$lvl[i] %in% c("domain","phylum")){
                                plotD.nodes$id[i]
                              }else if( plotD.nodes$id[i] %in% topGenus){
                                plotD.nodes$id[i]
                              }else NA
                            })


library(igraph)
net <- graph_from_data_frame(d=plotD.edges, vertices=plotD.nodes, directed=F) 
net

# subset the graph object 
sg1 <- decompose(net,mode="weak") 
sapply(sg1,length) 


# plot the bacteria sub network
net.sub <- sg1[[1]] 
plotName <- names(V(net.sub))[grep("d_",names(V(net.sub)))]


library("ggraph")
ggraph(net.sub, layout = 'igraph', algorithm="nicely") +
  geom_edge_link(aes(color=up.lvl)) +
  geom_node_point( aes(size = numARG, fill=lvl), shape=21) + #可以设定size： numARG, log10(numARG)
  geom_node_text(aes(label = labels), size=2, color="black", repel=T) +
  scale_edge_color_manual(values = taxlvl.cols[names(taxlvl.cols) %in% E(net.sub)$up.lvl]) +
  scale_fill_manual(values = taxlvl.cols) +
  scale_color_manual(values = c("white", "black")) +
  scale_size(range = c(3, 15)) +  # size = numARG时
  theme_void()

# this is a rough plot for the ARG host network, the refined layout are generated by Gephi and tailored on Adobe Illustrator.





