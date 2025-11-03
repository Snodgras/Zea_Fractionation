#!/usr/bin/env Rscript

library(tidyverse)

args<-commandArgs(T)
#args[1] = path to deletion bed file "exonic_M1_deletions.tsv"
#args[2] = "M1" or "M2"

print(paste("Reading in data files: ",Sys.time()))
#read in data file
del_data<-read_tsv(file = args[1],
                   col_names = c("CHROM_sb","Start_sb","Stop_sb","ID","QUAL","Strand","Genome","CHROM_del","Start_del","Stop_del","REF","ALT1","ALT2","Length_del","Genotype","Length_overlap"),
                   na = c("."))

#create extra variables
del_data<-mutate(del_data, Length_sb = Stop_sb - Start_sb)

if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                       args[2],"_del_reframe.tsv"))){

del_reframe<-del_data %>% group_by(ID, Genome) %>%
  reframe(ID = ID, Genome = Genome, Length_sb = Length_sb,
          total_overlap_del_bp=sum(Length_overlap),
          total_overlap_del_percent = (total_overlap_del_bp/Length_sb)*100) 

del_reframe<-select(del_reframe, ID, Length_sb) %>% unique() %>% 
  mutate(Length_sb_sizebin = case_when(Length_sb <= 10 ~ "1-10bp",
                                       Length_sb > 10 & Length_sb <= 100 ~ "11-100bp",
                                       Length_sb > 100 & Length_sb <= 1000 ~ "101-1000bp",
                                       Length_sb > 1000 ~ ">1000bp"),
         Length_sb_sizebin = factor(Length_sb_sizebin, levels = c("1-10bp","11-100bp","101-1000bp",">1000bp"))) %>%
  inner_join(y=del_reframe, by=c("ID","Length_sb","Genome"))

#write the del_reframe
print(paste("Writing deletion reframe tsv: ",Sys.time()))

write_tsv(del_reframe, paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                              args[2],"_del_reframe.tsv"))
}else{del_reframe <- read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                                     args[2],"_del_reframe.tsv"))}

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
                 "ZvTIL11" = "#018ce7",
                 "ZvTIL01" = "#018ce7",
                 "ZxTIL18" = "#036db2",
                 "ZxTIL25" = "#036db2",
                 "ZhRIMHU001" = "#2B5A78",
                 "ZdGigi" = "#433475",
                 "ZdMomo" = "#433475",
                 "ZnPI615697" = "#412151",
                 "TdFL" = "#F45B69",
                 "ZdGigi_4to1" = "#433475",
                 "ZdMomo_4to1" = "#433475",
                 "ZnPI615697_4to1" = "#412151",
                 "TdKS" = "#F45B69")

#Showing deletion overlap correlating with exon size:
print(paste("Plotting exon size by deleted segment length: ",Sys.time()))
if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
                       "_ExonLengthBins_by_PercentExonDelOverlap.png"))){
del_reframe %>% 
  filter(!is.na(Genome)) %>%
  group_by(Genome, Length_sb_sizebin) %>%
  summarize(mean.total_overlap_del_percent = mean(total_overlap_del_percent, na.rm = T),
            sd.total_overlap_del_percent = sd(total_overlap_del_percent, na.rm=T)) %>%
  inner_join(y=count(group_by(filter(del_reframe, !is.na(Genome)), Genome, Length_sb_sizebin))) %>%
  mutate(se.total_overlap_del_percent = sd.total_overlap_del_percent/sqrt(n),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                  "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                  "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                  "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                  "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                  "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" 
                                  ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
         Length_sb_sizebin = factor(Length_sb_sizebin, levels = c("1-10bp","11-100bp","101-1000bp",">1000bp")))%>%
  ggplot(aes(x=Length_sb_sizebin, y=mean.total_overlap_del_percent))+
  geom_bar(aes(fill=Genome), stat = "identity",position = "dodge")+
  geom_errorbar(aes(ymin=mean.total_overlap_del_percent-se.total_overlap_del_percent,
                    ymax=mean.total_overlap_del_percent+se.total_overlap_del_percent,
                    color=Genome), position = "dodge")+
  scale_fill_manual(values=genome_colors)+ scale_color_manual(values = rep("#000000",36))+
  guides(color="none")+
  theme_bw()+xlab("Exon length binned")+ylab("Mean Percent Exon Sequence Deleted")
ggsave(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
              "_ExonLengthBins_by_PercentExonDelOverlap.png"),
       device = "png",dpi=300,width = 8.5, height = 6)
}

#Do it by gene length:
print(paste("Plotting gene size by deleted segment length: ",Sys.time()))

if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
                       "_GeneLengthBins_by_PercentGeneDelOverlap.png"))){
dels_by_gene<-del_data %>% mutate(gene_ID = str_split(ID,";",simplify=T)[,2] %>% str_remove("Parent=")) %>%
  group_by(gene_ID) %>%
  summarize(gene_ID = gene_ID, Genome = Genome,Length_overlap=Length_overlap,
          Length_sb_gene = max(Stop_sb)-min(Start_sb)) %>%
  group_by(gene_ID, Genome) %>%
  mutate(total_overlap_del_bp=sum(Length_overlap),
         total_overlap_del_percent = (total_overlap_del_bp/Length_sb_gene)*100) %>% 
  select(-Length_overlap) %>% unique() %>% 
  mutate(gene_size = case_when(Length_sb_gene <= 1000 ~ "<= 1Kbp",
                               Length_sb_gene > 1000 & Length_sb_gene <= 2000 ~ "1-2Kbp",
                               Length_sb_gene > 2000 & Length_sb_gene <= 3000 ~ "2-3Kbp",
                               Length_sb_gene > 3000 & Length_sb_gene <= 4000 ~ "3-4Kbp",
                               Length_sb_gene > 4000 & Length_sb_gene <= 5000 ~ "4-5Kbp",
                               Length_sb_gene > 5000 ~ "> 5Kbp"),
         gene_size = factor(gene_size, levels = c("<= 1Kbp","1-2Kbp","2-3Kbp","3-4Kbp","4-5Kbp","> 5Kbp"))) %>% 
  filter(!is.na(Genome))

dels_by_gene %>% group_by(Genome, gene_size) %>%
  summarize(mean.total_overlap_del_percent = mean(total_overlap_del_percent, na.rm = T),
            sd.total_overlap_del_percent = sd(total_overlap_del_percent, na.rm=T)) %>%
  inner_join(y=count(group_by(dels_by_gene, Genome, gene_size))) %>%
  mutate(se.total_overlap_del_percent = sd.total_overlap_del_percent/sqrt(n),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                  "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                  "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                  "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                  "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                  "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" 
                                  ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
         gene_size = factor(gene_size, levels = c("1-2Kbp","2-3Kbp","3-4Kbp","4-5Kbp","> 5Kbp")))%>%
  ggplot(aes(x=gene_size, y=mean.total_overlap_del_percent))+
  geom_bar(aes(fill=Genome), stat = "identity",position = "dodge")+
  geom_errorbar(aes(ymin=mean.total_overlap_del_percent-se.total_overlap_del_percent,
                    ymax=mean.total_overlap_del_percent+se.total_overlap_del_percent,
                    color=Genome), position = "dodge")+
  scale_fill_manual(values=genome_colors)+ scale_color_manual(values = rep("#000000",36))+
  guides(color="none")+ theme_bw()+
  xlab("Gene length binned")+ylab("Mean Percent Gene Coding Sequence Deleted")
ggsave(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
              "_GeneLengthBins_by_PercentGeneDelOverlap.png"),
       device = "png",dpi=300,width = 8.5, height = 6)
}
#1:How many deletions overlap a single CDS?
filter(del_data, !is.na(Genome)) %>% group_by(ID, Genome) %>% count() %>% summary()
#read standard output

#2.How big are the overlaps?
print(paste("Plotting size of deleted segments: ",Sys.time()))
if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
                       "_DelOverlapSizes_barchart.png"))){
filter(del_data, !is.na(Genome)) %>%
  mutate(overlap_size_class = case_when(Length_overlap <= 10 ~ "1-10bp",
                                        Length_overlap > 10 & Length_overlap <=100 ~ "11-100bp",
                                        Length_overlap > 100 & Length_overlap <=1000 ~ "101-1000bp",
                                        Length_overlap > 1000~ ">1001bp"),
         overlap_size_class = factor(overlap_size_class, levels = c("1-10bp","11-100bp","101-1000bp",">1001bp"),
         Genome = factor(Genome,levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                 "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                  "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                  "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                  "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                  "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" 
                                  ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")))) %>%
  ggplot(aes(x=overlap_size_class, fill=Genome))+
  geom_bar(position = "dodge")+scale_fill_manual(values=genome_colors)+
  theme_bw()+xlab("Sizes of individual exonic deletions")+ylab("Counts")
ggsave(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
              "_DelOverlapSizes_barchart.png"),device="png",dpi=300,width = 8.5, height = 4)
}
#3. How many of the overlaps are multiples of 3?
print(paste("Plotting deleted segment multiple of 3 T/F: ",Sys.time()))
if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],"DelOverlapMult3_TFbars.png"))){
filter(del_data, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>%
  group_by(Genome, overlap_mult_3) %>% count() %>%
  mutate(percent = (n/69269)*100) %>%
  ggplot(aes(x=Genome, y=percent, fill=overlap_mult_3))+
  geom_bar(position="fill",stat = "identity")+
  theme_bw()+xlab("")+ylab("Proportion")+scale_fill_manual(name="Del. exon sequence \nlength is multiple of 3bp",
                                                           values=c("grey","black"))
ggsave(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],"DelOverlapMult3_TFbars.png")
       ,device="png",dpi=300,width = 8.5, height = 8.5)
}
#3.a How many CDS only have multiple of 3 deletions vs. mix vs. just non-mult. of 3 deletions?
if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                       args[2],"_mult_3_deletions.tsv"))){
mult_3_dels<-filter(del_data, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>% 
  select(ID, Genome, overlap_mult_3) %>% unique() %>% 
  group_by(ID, Genome) %>%
  reframe(ID = ID, Genome = Genome, 
          Deletion_Types = case_when(all(isTRUE(overlap_mult_3)) ~ "Only mult 3",
                                     all(isFALSE(overlap_mult_3)) ~ "Only non-mult 3",
                                     .default = "Mixed"))
print(paste("Writing multiple of 3 deletion tsv: ",Sys.time()))

write_tsv(mult_3_dels,paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                 args[2],"_mult_3_deletions.tsv"))
}else{mult_3_dels<-read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                                   args[2],"_mult_3_deletions.tsv"))}
print(paste("Plotting multiple of 3 deletion types by genome: ",Sys.time()))
if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
                       "_Mult3DelTypes_barchart.png"))){
mult_3_dels %>%
mutate(Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                          "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                          "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" 
                                          ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
  Deletion_Types = factor(Deletion_Types, levels=c("Only mult 3","Mixed","Only non-mult 3"))) %>%
  ggplot(aes(x=Deletion_Types, fill=Genome))+
  geom_bar(stat = "count", position="dodge")+
  xlab("Deletion Types")+theme_bw()+
  scale_fill_manual(values=genome_colors)
ggsave(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Plots/",args[2],
              "_Mult3DelTypes_barchart.png"),device="png",dpi=300,width = 8.5, height = 4)
}

if(!file.exists(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                       args[2],"_mult_3_deletions.genes.tsv"))){
filter(del_data, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>% 
  select(ID, Genome, overlap_mult_3) %>% unique() %>% 
  group_by(ID, Genome) %>%
    mutate(Deletion_Types = if(all(isTRUE(overlap_mult_3))){"Only mult 3"}else{
      if(all(isFALSE(overlap_mult_3))){"Only non-mult 3"}else{
        "Mixed"}}) %>%
    select(-overlap_mult_3) %>% 
    unique() %>% 
  mutate(gene_ID = str_split(ID, ";", simplify=TRUE)[,2] %>% str_remove("Parent=")) %>%
  group_by(gene_ID, Genome) %>% 
  mutate(gene_Deletion_Types = if(length(unique(Deletion_Types)) > 1){"Mixed"}else{unique(Deletion_Types)}) %>%
  write_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/",
                   args[2],"_mult_3_deletions.genes.tsv"))
}