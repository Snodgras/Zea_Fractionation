library(ggplot2)
library(compiler)
library(tidyverse)
enableJIT(3)
#library(ggplot2)
#library("Cairo")
changetoM <- function ( position ){
  position=position/1000000; 
  paste(position, "M", sep="")
}
#For the perl script
#data=read.table("Zv-TIL01.Sbicolor.tab")
#data = data[which(data$V1 %in% c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")),]
#data = data[which(data$V3 %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10")),]
#data$V1 = factor(data$V1, levels=c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"))
#data$V3 = factor(data$V3, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"))
#ggplot(data = data, aes(x=V4, y=V2))+
#  geom_point(size=0.5, aes(color=V5))+
#  facet_grid(V1~V3, scales="free", space="free") + 
##  theme_grey(base_size=120)+
#  labs(x="Zv-TIL01", y="Sbicolor_313")+
#  scale_x_continuous(labels=changetoM)+
#  scale_y_continuous(labels=changetoM)+
#  theme(axis.line = element_blank(),
#        panel.background = element_blank(),
#        panel.border = element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
#        axis.text.y = element_text(colour = "black"),
#        legend.position="none",
#        axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
#ggsave("Zv-TIL01.Sbicolor_313.dotplot.pdf",device = "pdf")

#For the anchors
data=read.table("Anchors/Sb313_TdFL.bin2_anchorwave.anchors",head=TRUE)
data = data[which(data$refChr %in% c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")),]
#for split genomes:
data$queryChr<-str_split(string = data$queryChr, pattern = ":", simplify = TRUE)[,1]
data = data[which(data$queryChr %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18")),]
data$refChr = factor(data$refChr, levels=c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"))
data$queryChr = factor(data$queryChr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
ggplot(data = data, aes(x=queryStart, y=referenceStart))+
  geom_point(size=0.5, aes(color=strand))+
  facet_grid(refChr~queryChr, scales="free", space="free") + 
#  theme_grey(base_size=120)+
  labs(x="Av", y="Sbicolor_313")+
  scale_x_continuous(labels=changetoM)+
  scale_y_continuous(labels=changetoM)+
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
       panel.border = element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
        axis.text.y = element_text(colour = "black"),
       legend.position="none",
        axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
#ggsave("Av.Sb313.anchor.dotplot.pdf",device = "pdf")

makeDotplot<-function(df, GenomeName, Split){
  data = read.table(df, header = TRUE)
  data = data[which(data$refChr %in% c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")),]
  if(Split == TRUE){
    data$queryChr<-str_split(string = data$queryChr, pattern = ":", simplify = TRUE)[,1]
  }
  data = data[which(data$queryChr %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18")),]
  data$refChr = factor(data$refChr, levels=c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"))
  data$queryChr = factor(data$queryChr, levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"))
  Dotplot<-ggplot(data = data, aes(x=queryStart, y=referenceStart))+
    geom_point(size=0.5, aes(color=strand))+
    facet_grid(refChr~queryChr, scales="free", space="free") + 
    #  theme_grey(base_size=120)+
    labs(x=GenomeName, y="Sbicolor_313")+
    scale_x_continuous(labels=changetoM)+
    scale_y_continuous(labels=changetoM)+
    theme(axis.line = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
          axis.text.y = element_text(colour = "black"),
          legend.position="none",
          axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
  return(Dotplot)
}
makeDotplot("Anchors/Sb313_Av_anchorwave.anchors", "Av", FALSE)
ggsave("Sb313vsAv.anchor.dotplot.pdf", device="pdf")
#go to Anchors director and do: 
#for i in *.anchors ; do echo Anchors/$i >> filelist.txt ; done

file.list<-read.delim("Anchors/filelist.txt",header = FALSE)
str(file.list)
file.list<-mutate(file.list, genome = str_split(string = V1, pattern = "_", simplify = TRUE)[,2])

for(i in 1:nrow(file.list)){
  makeDotplot(file.list$V1[i], file.list$genome[i], TRUE)
  ggsave(paste0("Sb313vs",file.list$genome[i],".anchor.dotplot.pdf"), device = "pdf")
}

anchor.counts<-read_tsv("Anchors/anchor.counts",col_names = c("Genome","Count_Sobic","Count_localAlignment"))
anchor.counts<-mutate(anchor.counts, Count_All = Count_Sobic +Count_localAlignment,
                      Percent_Sobic = Count_Sobic/34211,
                      Percent_localAlignment = Count_localAlignment/34211,
                      Percent_All = Count_All/34211)
ggplot(anchor.counts, aes(x=Genome, y=Percent_All))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

anchor.counts$Genome<-anchor.counts$Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                                               "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                                               "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                                               "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZlRIL003","TdFL", "Av"))
genome_colors<-c("ZmB73" = "#DAFDD8",
                 "ZmB97" = "#DAFDD8",
                 "ZmCML103" = "#7CCBB2",
                 "ZmCML228" = "#7CCBB2",
                 "ZmCML247" = "#7CCBB2",
                 "ZmCML277" = "#7CCBB2",
                 "ZmCML322" = "#7CCBB2",
                 "ZmCML333" = "#7CCBB2",
                 "ZmCML52" = "#7CCBB2",
                 "ZmCML69" = "#7CCBB2",
                 "ZmHP301" = "#92DDB0",
                 "ZmIL14H" = "#DAFDD8",
                 "ZmKi11" = "#7CCBB2",
                 "ZmKi3" = "#7CCBB2",
                 "ZmKy21" = "#DAFDD8",
                 "ZmM162W" = "#DAFDD8",
                 "ZmM37W" = "#92DDB0",
                 "ZmMo18W" = "#92DDB0",
                 "ZmMS71" = "#DAFDD8",
                 "ZmNC350" = "#7CCBB2",
                 "ZmNC358" = "#7CCBB2",
                 "ZmOh43" = "#DAFDD8",
                 "ZmOh7b" = "#DAFDD8",
                 "ZmP39" = "#DAFDD8",
                 "ZmTx303" = "#92DDB0",
                 "ZmTzi8" = "#7CCBB2",
                 "ZvTIL11" = "#4DA89D",
                 "ZvTIL01" = "#4DA89D",
                 "ZxTIL18" = "#03A0B5",
                 "ZxTIL25" = "#03A0B5",
                 "ZhRIMHU001" = "#247590",
                 "ZdGigi" = "#2B5A78",
                 "ZdMomo" = "#2B5A78",
                 "ZnPI615697" = "#114B5F",
                 "TdFL" = "#F45B69",
                 "Av" = "black",
                 "ZlRIM003"="gray")
ggplot(anchor.counts, aes(x=Genome, y=Percent_Sobic))+
  #geom_point(aes(color=Genome))+
  geom_point()+
  theme_bw()+
  #scale_color_manual(values=genome_colors)+
  theme(axis.text.x = element_text(angle=90))

