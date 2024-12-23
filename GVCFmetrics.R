#QC GVCF Alignments with the VCF Metrics results from Tassel
library(tidyverse)

#in bash do: grep "ALL" GVCFmetrics.tsv > clean.GVCFmetrics.tsv
#grep "ALL" GVCFmetrics.with4to1andTdKS.tsv > clean.GVCFmetrics.with4to1andTdKS.tsv

GVCF_metrics<-read_tsv("clean.GVCFmetrics.tsv", 
                       col_names = c("name","numSNPs","numIndels","numIns",
                                     "numDel","numBasesInserted","numBasesDeleted","numBasesInGaps",
                                     "totalBases","percentRefAligned","percentRefMapped","percentRefDeleted",
                                     "percentRefInserted","meanIndelSize","lowerQuartileIndelSize","medianIndelSize",
                                     "upperQuartileIndelSize","largestIndel","meanInsSize","lowerQuartileInsSize",
                                     "medianInsSize","upperQuartileInsSize","largestInsSize","meanDelSize",
                                     "lowerQuartileDelSize","medianDelSize","upperQuartileDelSize","largestDelSize",
                                     "extra"))
GVCF4to1andTdKS_metrics<-read_tsv("clean.GVCFmetrics.with4to1andTdKS.tsv", 
                                  col_names = c("name","numSNPs","numIndels","numIns",
                                                "numDel","numBasesInserted","numBasesDeleted","numBasesInGaps",
                                                "totalBases","percentRefAligned","percentRefMapped","percentRefDeleted",
                                                "percentRefInserted","meanIndelSize","lowerQuartileIndelSize","medianIndelSize",
                                                "upperQuartileIndelSize","largestIndel","meanInsSize","lowerQuartileInsSize",
                                                "medianInsSize","upperQuartileInsSize","largestInsSize","meanDelSize",
                                                "lowerQuartileDelSize","medianDelSize","upperQuartileDelSize","largestDelSize",
                                                "extra"))

#get rid of column "extra" and split the name by genome and refChr.bin
clean_metrics<-GVCF_metrics %>% select(-extra) %>% 
  mutate(genome = str_split(string=name, pattern="_",simplify=T)[,2],
         refChr.bin = str_split(string=name, pattern="_",simplify=T)[,3],
         refChr = str_split(string = refChr.bin, pattern="[.]",simplify = T)[,1],
         bin = str_split(string = refChr.bin,pattern="[.]",simplify = T)[,2]) 

genome_colors<-c("ZmB73" = "#92ddb0",
                 "ZmB97" = "#92ddb0",
                 "ZmCML103" = "#4da89d",
                 "ZmCML228" = "#4da89d",
                 "ZmCML247" = "#4da89d",
                 "ZmCML277" = "#4da89d",
                 "ZmCML322" = "#4da89d",
                 "ZmCML333" = "#4da89d",
                 "ZmCML52" = "#4da89d",
                 "ZmCML69" = "#4da89d",
                 "ZmHP301" = "#92ddb0",
                 "ZmIL14H" = "#92ddb0",
                 "ZmKi11" = "#4da89d",
                 "ZmKi3" = "#4da89d",
                 "ZmKy21" = "#92ddb0",
                 "ZmM162W" = "#92ddb0",
                 "ZmM37W" = "#7ccbb2",
                 "ZmMo18W" = "#7ccbb2",
                 "ZmMS71" = "#92ddb0",
                 "ZmNC350" = "#4da89d",
                 "ZmNC358" = "#4da89d",
                 "ZmOh43" = "#92ddb0",
                 "ZmOh7b" = "#92ddb0",
                 "ZmP39" = "#92ddb0",
                 "ZmTx303" = "#7ccbb2",
                 "ZmTzi8" = "#4da89d",
                 "ZvTIL11" = "#03a0b5",
                 "ZvTIL01" = "#03a0b5",
                 "ZxTIL18" = "#247590",
                 "ZxTIL25" = "#247590",
                 "ZhRIMHU001" = "#2B5A78",
                 "ZdGigi" = "#433475",
                 "ZdMomo" = "#433475",
                 "ZnPI615697" = "#412151",
                 "TdFL" = "#F45B69")
chr_colors<-c("Chr01"="#5E0D59", "Chr02"="#8A1465", "Chr03"="#B61B5E", 
              "Chr04"="#E12346", "Chr05"="#ED524A", "Chr06"="#EE9E7C",
              "Chr07"="#FFB185","Chr08"="#FC9E40","Chr09"="#FEB820",
              "Chr10"="#FFD53D")
clean_metrics$genome<-clean_metrics$genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                              "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                              "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                              "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","TdFL"))
genome_colors4to1<-c("ZmB73" = "#92ddb0", "ZmB97" = "#92ddb0","ZmCML103" = "#4da89d","ZmCML228" = "#4da89d","ZmCML247" = "#4da89d","ZmCML277" = "#4da89d","ZmCML322" = "#4da89d","ZmCML333" = "#4da89d","ZmCML52" = "#4da89d","ZmCML69" = "#4da89d","ZmHP301" = "#92ddb0", "ZmIL14H" = "#92ddb0","ZmKi11" = "#4da89d",
                 "ZmKi3" = "#4da89d","ZmKy21" = "#92ddb0","ZmM162W" = "#92ddb0","ZmM37W" = "#7ccbb2","ZmMo18W" = "#7ccbb2","ZmMS71" = "#92ddb0","ZmNC350" = "#4da89d","ZmNC358" = "#4da89d","ZmOh43" = "#92ddb0","ZmOh7b" = "#92ddb0","ZmP39" = "#92ddb0","ZmTx303" = "#7ccbb2",
                 "ZmTzi8" = "#4da89d","ZvTIL11" = "#03a0b5","ZvTIL01" = "#03a0b5","ZxTIL18" = "#247590","ZxTIL25" = "#247590","ZhRIMHU001" = "#2B5A78",
                 "ZdGigi" = "#433475",
                 "ZdMomo" = "#433475",
                 "ZnPI615697" = "#412151",
                 "TdFL" = "#F45B69",
                 "ZdGigi.4to1" = "#433475",
                 "ZdMomo.4to1" = "#433475",
                 "ZnPI615697.4to1" = "#412151",
                 "TdKS" = "#F45B69")
genome_shape<-c("ZmB73" = 1, "ZmB97" = 1,"ZmCML103" = 1,"ZmCML228" = 1,"ZmCML247" = 1,"ZmCML277" = 1,"ZmCML322" = 1,"ZmCML333" = 1,"ZmCML52" = 1,"ZmCML69" = 1,"ZmHP301" = 1, "ZmIL14H" = 1,"ZmKi11" = 1,
                "ZmKi3" = 1,"ZmKy21" = 1,"ZmM162W" = 1,"ZmM37W" = 1,"ZmMo18W" = 1,"ZmMS71" = 1,"ZmNC350" = 1,"ZmNC358" = 1,"ZmOh43" = 1,"ZmOh7b" = 1,"ZmP39" = 1,"ZmTx303" = 1,
                "ZmTzi8" = 1,"ZvTIL11" = 1,"ZvTIL01" = 1,"ZxTIL18" = 1,"ZxTIL25" = 1,"ZhRIMHU001" = 1,
                "ZdGigi" = 1,
                "ZdMomo" = 1,
                "ZnPI615697" = 1,
                "TdFL" = 1,
                "ZdGigi.4to1" = 17,
                "ZdMomo.4to1" = 17,
                "ZnPI615697.4to1" = 17,
                "TdKS" = 17)
GVCF4to1andTdKS_metrics$name<-str_replace_all(string=GVCF4to1andTdKS_metrics$name, pattern = "_4to1", replacement = ".4to1")
clean4to1andTdKS_metrics<-GVCF4to1andTdKS_metrics %>% select(-extra) %>% 
  mutate(genome = str_split(string=name, pattern="_",simplify=T)[,2],
         refChr.bin = str_split(string=name, pattern="_",simplify=T)[,3],
         refChr = str_split(string = refChr.bin, pattern="[.]",simplify = T)[,1],
         bin = str_split(string = refChr.bin,pattern="[.]",simplify = T)[,2]) 

clean4to1andTdKS_metrics$genome<-clean4to1andTdKS_metrics$genome %>% 
  factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,
                  "ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                  "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                  "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                  "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,
                  "ZmNC358" ,"ZmTzi8" ,"ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,
                  "ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi.4to1","ZdMomo.4to1","ZnPI615697.4to1","TdFL","TdKS"))

ggplot(clean4to1andTdKS_metrics, aes(color=genome, shape = genome, x=bin, y=numSNPs))+
  geom_jitter(height = 0)+
  facet_wrap(vars(refChr))+
  theme_bw()+
  scale_color_manual(values=genome_colors4to1)+scale_shape_manual(values = genome_shape)
filter(clean4to1andTdKS_metrics, genome %in% c("ZdGigi","ZdMomo","ZnPI615697" ,
                "TdFL","ZdGigi.4to1","ZdMomo.4to1","ZnPI615697.4to1","TdKS")) %>%
  ggplot(aes(color=genome, shape = genome, x=bin, y=numSNPs))+
  geom_jitter(height = 0)+
  facet_wrap(vars(refChr))+
  theme_bw()+
  scale_color_manual(values=genome_colors4to1)+scale_shape_manual(values = genome_shape)


plot_4to1_metrics<-function(metric, title){
  p<-filter(clean4to1andTdKS_metrics, genome %in% c("ZdGigi","ZdMomo","ZnPI615697" ,
                                                   "TdFL","ZdGigi.4to1","ZdMomo.4to1","ZnPI615697.4to1","TdKS")) %>%
    ggplot(aes_string(color="genome", shape = "genome", x="bin", y=metric))+
    geom_jitter(height = 0)+
    facet_wrap(vars(refChr))+
    theme_bw()+
    scale_color_manual(values=genome_colors4to1)+scale_shape_manual(values = genome_shape)+
    ggtitle(title)+
    ylab("")+xlab("")
  return(p)
}

plot_4to1_metrics("numSNPs","Number of SNPs in GVCF")
ggsave("GVCFmetrics.4to1.numSNPs.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("numIns","Number of Insertions in GVCF")
ggsave("GVCFmetrics.4to1.numIns.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("numDel","Number of Deletions in GVCF")
ggsave("GVCFmetrics.4to1.numDel.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("percentRefAligned","Percent Ref Aligned in GVCF")
ggsave("GVCFmetrics.4to1.percentRefAligned.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("percentRefMapped","Percent Ref Mapped in GVCF")
ggsave("GVCFmetrics.4to1.percentRefMapped.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("percentRefDeleted","Percent Ref Deleted in GVCF")
ggsave("GVCFmetrics.4to1.percentRefDeleted.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("meanIndelSize","Mean Indel Size in GVCF")
ggsave("GVCFmetrics.4to1.meanIndelSize.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("percentRefInserted","Percent Ref Inserted in GVCF")
ggsave("GVCFmetrics.4to1.percentRefInserted.pdf",device="pdf",dpi=300, width = 5.5)

plot_4to1_metrics("largestIndel","Largest Indel in GVCF")
ggsave("GVCFmetrics.4to1.largestIndel.pdf",device="pdf",dpi=300, width = 5.5)

#Plotting the 4 to 1 AW runs and TdKS along with the 2 to 1 runs

ggplot(clean4to1andTdKS_metrics, aes(color=genome, x=bin, y=numSNPs))+
  geom_jitter(height = 0)+
  facet_wrap(vars(refChr))+
  theme_bw()+
  scale_color_manual(values=genome_colors)+
  ggtitle("Number of SNPs")+
  ylab("")+xlab("")

plot_metrics<-function(metric, title){
  p<-ggplot(clean4to1andTdKS_metrics,aes_string(color="genome", shape = "genome", x="bin", y=metric))+
    geom_jitter(height = 0)+
    facet_wrap(vars(refChr))+
    theme_bw()+
    scale_color_manual(values=genome_colors4to1)+scale_shape_manual(values = genome_shape)+
    ggtitle(title)+
    ylab("")+xlab("")
  return(p)
}
plot_metrics("numSNPs","Number of SNPs in GVCF")
ggsave("GVCFmetrics.4to1.numSNPs.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("numIns","Number of Insertions in GVCF")
ggsave("GVCFmetrics.numIns.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("numDel","Number of Deletions in GVCF")
ggsave("GVCFmetrics.numDel.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("percentRefAligned","Percent Ref Aligned in GVCF")
ggsave("GVCFmetrics.percentRefAligned.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("percentRefMapped","Percent Ref Mapped in GVCF")
ggsave("GVCFmetrics.percentRefMapped.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("percentRefDeleted","Percent Ref Deleted in GVCF")
ggsave("GVCFmetrics.percentRefDeleted.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("meanIndelSize","Mean Indel Size in GVCF")
ggsave("GVCFmetrics.meanIndelSize.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("percentRefInserted","Percent Ref Inserted in GVCF")
ggsave("GVCFmetrics.percentRefInserted.pdf",device="pdf",dpi=300, width = 5.5)

plot_metrics("largestIndel","Largest Indel in GVCF")
ggsave("GVCFmetrics.largestIndel.pdf",device="pdf",dpi=300, width = 5.5)

####Plotting 4 to 1 metrics without the split by Ref Chromosome
filter(clean4to1andTdKS_metrics, genome %in% c("ZdGigi","ZdMomo","ZnPI615697" ,
                                               "TdFL","ZdGigi.4to1","ZdMomo.4to1","ZnPI615697.4to1","TdKS")) %>%
  ggplot(aes(color=genome,
             #fill = genome,
             shape = genome, 
             x=bin, y=percentRefMapped))+
  geom_jitter(height = 0, width = 0.25)+
  facet_wrap(vars(refChr))+
  #geom_boxplot()+
  theme_bw()+
  scale_color_manual(values = c("TdKS" = "darkred", "TdFL"="red",
                                "ZnPI615697.4to1"="darkorchid4", "ZnPI615697"="mediumorchid",
                                "ZdMomo.4to1"="royalblue4", "ZdMomo"="royalblue",
                                "ZdGigi.4to1"="springgreen4","ZdGigi"="seagreen3"))+scale_shape_manual(values = genome_shape)+
  ggtitle("4 to 1 Mapping + TdKS % Ref Mapped")+
  ylab("")+xlab("")

library(ggridges)
ggplot(clean4to1andTdKS_metrics, aes(x=percentRefMapped, y=genome, fill = genome))+
  geom_density_ridges()+
  facet_wrap(~bin)+
  scale_fill_manual(values = c("TdKS" = "deeppink4", "TdFL"="hotpink",
                               "ZnPI615697.4to1"="darkorchid4", "ZnPI615697"="mediumorchid",
                               "ZdMomo.4to1"="royalblue4", "ZdMomo"="royalblue",
                               "ZdGigi.4to1"="springgreen4","ZdGigi"="seagreen3",
                               "ZmB73" = "darkgrey", "ZmB97" = "darkgrey","ZmCML103" = "darkgrey","ZmCML228" = "darkgrey","ZmCML247" = "darkgrey","ZmCML277" = "darkgrey","ZmCML322" = "darkgrey","ZmCML333" = "darkgrey","ZmCML52" = "darkgrey","ZmCML69" = "darkgrey","ZmHP301" = "darkgrey", "ZmIL14H" = "darkgrey","ZmKi11" = "darkgrey",
                               "ZmKi3" = "darkgrey","ZmKy21" = "darkgrey","ZmM162W" = "darkgrey","ZmM37W" = "darkgrey","ZmMo18W" = "darkgrey","ZmMS71" = "darkgrey","ZmNC350" = "darkgrey","ZmNC358" = "darkgrey","ZmOh43" = "darkgrey","ZmOh7b" = "darkgrey","ZmP39" = "darkgrey","ZmTx303" = "darkgrey",
                               "ZmTzi8" = "darkgrey","ZvTIL11" = "darkgrey","ZvTIL01" = "darkgrey","ZxTIL18" = "darkgrey","ZxTIL25" = "darkgrey","ZhRIMHU001" = "darkgrey"))
