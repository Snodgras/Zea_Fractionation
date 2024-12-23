#This is to build the phylogeny figure from Orthofinder Species Tree
#install.packages("ape")
library(BiocManager)
BiocManager::install("ggtree")
library(ape)
library(ggtree)
library(tidyverse)

tree <- read.tree(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/orthofinder/Results_Aug21/Species_Tree/SpeciesTree_rooted.txt")
treewNodes<-read.tree(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/orthofinder/Results_Aug21/Species_Tree/SpeciesTree_rooted_node_labels.txt")
#sppcode <- read.delim("/work/mash-covid/genespace/temp/sppcode.tsv")
#tip <- tree$tip.label
#myLabels <- sppcode[sppcode$sixLetterCode %in% tip ,]
#myLabels <- myLabels %>% select(sixLetterCode, genusSpp)
myLabels<-tibble(genussSpp = c("Anatherum virginicum","Sorghum bicolor","Tripsacum dactyloides","Zea nicaraguensis","Zea diploperennis","Zea diploperennis",
                               "Zea huehuetenangensis","Zea mays mexicana","Zea mays mexicana","Zea mays parviglumis","Zea mays parviglumis",
                               "Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays",
                               "Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays",
                               "Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays",
                               "Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays","Zea mays mays"
                               ),
                 label = tree$tip.label)
myLabels
p <- ggtree(tree)
#p  %<+% myLabels  +  geom_tiplab(
p + geom_tiplab(
  align = F,
  size = 3,
  linesize = .3,
 # aes(fill = genusSpp),
  color = "black",
  geom = "label", fontface='italic',
  label.padding = unit(0.15, "lines"),
  label.size = 0) +       theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.key = element_blank()
  ) + geom_treescale(x = 0.2, y = 10) 
+ geom_nodelab(
    size = 2,
    color = "black",
    geom = "label", fill = "lightgreen", fontface='bold')
ggsave("diploids_tree.png", width = 10, height = 5)

ggtree(treewNodes)+geom_tiplab(
  align = F,
  size = 3,
  linesize = .3,
  # aes(fill = genusSpp),
  color = "black",
  geom = "label", fontface='italic',
  label.padding = unit(0.15, "lines"),
  label.size = 0) +       theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.key = element_blank()
  ) + geom_treescale(x = 0.2, y = 10) + geom_nodelab(
  size = 2,
  color = "black",
  geom = "label", fontface='bold', fill = "lightgreen")

