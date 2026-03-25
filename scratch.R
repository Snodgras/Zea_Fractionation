#Scratch

## Double checking that fractionation calls aren't completely correlated with strand
library(tidyverse)

full.fractionation.status %>% colnames()
ref_Sb313.cds<-read_tsv("ref_Sb313.cds.bed", col_names = c("Chr","Start","Stop","ID", "QUAL","Strand"))

test<-left_join(full.fractionation.status, ref_Sb313.cds, by = "ID") %>% 
  select(Strand, ends_with("M1"),ends_with("M2")) %>% 
  mutate(Strand = case_when(Strand == "-" ~ 1, Strand == "+" ~ 0))

placeholder<-select(test, Strand, TdFL.M1) %>% na.omit()
t<-cor.test(placeholder$Strand,placeholder$TdFL.M1)

colnames(test)[2:79]

results<-tibble(cor_estimate = NA, 
                Genome_Subgenome = NA, 
                p.value = NA,
                lower.CI = NA,
                upper.CI = NA)
for(g in 2:79){
  placeholder<-test[,c(1,g)] %>% na.omit()
  t<-cor.test(placeholder$Strand,pull(placeholder[,2]))
  results<-add_row(results, cor_estimate = t$estimate, 
                      Genome_Subgenome = colnames(test)[g], 
                      p.value = t$p.value,
                      lower.CI = t$conf.int[1],
                      upper.CI = t$conf.int[2])
}
results<-na.omit(results)
filter(results, p.value <= 0.05) %>% nrow() #19 / 78


full.fractionation.status %>% 
  filter(str_detect(ID, "Sobic.010G")) %>% 
  select(ZmB73.M1) %>% group_by(ZmB73.M1) %>% count()

ZmB73_M1_10_crossmap<-read_tsv("ZmB73_Chr10.bin1.crossmap.tsv",
                               col_names = c("Chr_Zm","Start_Zm","Stop_Zm","ID","QUAL","Strand", "Map_ratio"))
ZmB73_M1_10_unmap<-read_tsv("ZmB73_Chr10.bin1.crossmap.tsv.unmap",
                            col_names = c("Chr_Sb","Start_Sb","Stop_Sb","ID","QUAL","Strand","Fail", "Map_ratio"))
Chr10_IDs<-full.fractionation.status %>% 
  filter(str_detect(ID, "Sobic.010G")) %>% select(ID) %>% pull()

comparison<-tibble(OG_call = NA, CM_call = NA, ID=NA)
for(i in Chr10_IDs){
  if(i %in% ZmB73_M1_10_crossmap$ID){
    comparison<-add_row(comparison,
                        OG_call = filter(full.fractionation.status, ID == i) %>% select(ZmB73.M1) %>% pull(),
                        CM_call = filter(ZmB73_M1_10_crossmap, ID == i) %>% select(Map_ratio) %>% pull(),
                        ID = i)
  }else{
    if(i %in% ZmB73_M1_10_unmap$ID){
      comparison<-add_row(comparison,
                          OG_call = filter(full.fractionation.status, ID == i) %>% select(ZmB73.M1) %>% pull(),
                          CM_call = filter(ZmB73_M1_10_unmap, ID == i) %>% select(Map_ratio) %>% pull(),
                          ID =i)
      
    }
    else{
      comparison<-add_row(comparison,
                          OG_call = filter(full.fractionation.status, ID == i) %>% select(ZmB73.M1) %>% pull(),
                          CM_call = NA,
                          ID = i)
      
    }
    
  }
  
}
comparison<-comparison[-1,]
comparison<-mutate(comparison, CM_call = str_remove_all(CM_call, "map_ratio="))
comparison<-comparison %>% 
  mutate(CM_summary = case_when(#CM_call == "Unmap" ~ NA,
                                CM_call == "Unmap" ~ 1,
                                CM_call >= 1.0000 ~ 0,
                                CM_call < 1 ~ 1)) 
ggplot(comparison,aes(x=OG_call, y=CM_summary)) +
  geom_jitter()+
  theme_bw()

filter(comparison, OG_call != CM_summary) %>% View()
25/nrow(comparison)         
filter(comparison, CM_call == "Unmap") %>% group_by(OG_call) %>% count()

#Unmap is if there's a threshhold 

#Try to figure out what's causing the disagreement 
questionable_ids<-filter(comparison, OG_call != CM_summary) %>% select(ID) %>% pull()

filter(ZmB73_M1_10_crossmap, ID %in% questionable_ids) %>% view()

comparison %>% 
  filter(OG_call == CM_summary) %>%
  select(ID) %>% 
  pull() %>% write_lines("test_B73_SbChr10_M1.cds.ids.txt")
#grep -f test_B73_SbChr10_M1.cds.ids.txt ZmB73_Chr10.bin1.crossmap.tsv > test_B73_SbCh10_M1.cds.B73coords.bed

#deletion bed files test
del_test<-read_tsv("Del_test/test.tsv", 
                   col_names = c("CHROM_sb","Start_sb","Stop_sb","ID","QUAL","Strand","Genome","CHROM_del","Start_del","Stop_del","REF","ALT1","ALT2","Length_del","Genotype","Length_overlap"),
                   na = c("."))

del_test<-mutate(del_test, Length_sb = Stop_sb - Start_sb) 

#this plot shows that by using the total overlap bps, you don't have more deleted bps than there are basepairs in the exon
del_test %>% group_by(ID, Genome) %>%
  reframe(ID = ID, Genome = Genome, Length_sb = Length_sb,
            total_overlap_del_bp=sum(Length_overlap)) %>%
  ggplot(aes(x= Length_sb, y=total_overlap_del_bp))+
  geom_point(aes(color = Genome), alpha = 0.5)+
  geom_abline(slope=1, intercept = 0)+
  theme_bw()+
  xlab("Exon length (bp)")+ylab("Total bp overlapping with deletions")

del_test_reframe<-del_test %>% group_by(ID, Genome) %>%
  reframe(ID = ID, Genome = Genome, Length_sb = Length_sb,
          total_overlap_del_bp=sum(Length_overlap),
          total_overlap_del_percent = (total_overlap_del_bp/Length_sb)*100) 

#what's the average exon length in each decile?
cds_size<-select(del_test_reframe, ID, Length_sb) %>% unique() %>% 
  mutate(Length_sb_sizebin = case_when(Length_sb <= 10 ~ "1-10bp",
                                       Length_sb > 10 & Length_sb <= 100 ~ "11-100bp",
                                       Length_sb > 100 & Length_sb <= 1000 ~ "101-1000bp",
                                       Length_sb > 1000 ~ ">1000bp"),
         Length_sb_sizebin = factor(Length_sb_sizebin, levels = c("1-10bp","11-100bp","101-1000bp",">1000bp")))
ggplot(cds_size, aes(x=Length_sb_sizebin, y=Length_sb))+
  geom_boxplot()+
  theme_bw()+xlab("Size bins of Exon Length")+ylab("Length of Exon (bp)")
ggsave("Del_test/ExonLengthSizeBins_by_ExonLengthBP.png", device = "png",dpi=300)

del_test_reframe<-inner_join(x=del_test_reframe,y=cds_size, by=c("ID","Length_sb"))

#testing out various plots for utility
del_test_reframe %>% 
  filter(!is.na(Genome)) #%>%
  #ggplot(aes(x=Length_sb_decile, y=total_overlap_del_percent))+
  #geom_jitter(height=0, aes(color=Genome))+
  #geom_boxplot(aes(fill = Genome))+
  #ggplot(aes(x=Length_sb, y=total_overlap_del_percent))+
  #geom_point(aes(color=Genome))+
  #ggplot(aes(x=total_overlap_del_percent, after_stat(count)))+
  #geom_histogram(binwidth = 1)+facet_wrap(vars(Genome))+
  #geom_density(aes(fill=Genome), alpha=0.5)+
  #theme_bw()

#Showing deletion overlap correlating with exon size:
del_test_reframe %>% 
  filter(!is.na(Genome)) %>%
  ggplot(aes(x=Length_sb_sizebin, y=total_overlap_del_percent))+
  geom_boxplot(aes(fill = Genome))+
  theme_bw()+xlab("Size bins of exon length")+ylab("Percent of exon overlapping with deletions")
ggsave("Del_test/ExonLengthDeciles_by_PercentExonDelOverlap.png",device = "png",dpi=300)

#Do it by gene length:
del_test_gene<-del_test %>% mutate(gene_ID = str_split(ID,";",simplify=T)[,2] %>% str_remove("Parent=")) %>%
  group_by(gene_ID) %>%
  reframe(gene_ID = gene_ID, Genome = Genome,Length_overlap=Length_overlap,
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
         gene_size = factor(gene_size, levels = c("<= 1Kbp","1-2Kbp","2-3Kbp","3-4Kbp","4-5Kbp","> 5Kbp")))
del_test_gene %>% filter(!is.na(Genome)) %>% 
  ggplot(aes(x=gene_size, y=total_overlap_del_percent,fill=Genome))+
  geom_boxplot()+
  theme_bw()


#1. How many deletions overlap a single CDS?
#2. How big are the overlaps? 
#3. How many of the overlaps are multiples of 3?

#1:
filter(del_test, !is.na(Genome)) %>% group_by(ID, Genome) %>% count() %>%
  ggplot(aes(x=Genome, y=n))+
  geom_boxplot(aes(color=Genome))+
  theme_bw()+xlab("")+ylab("Number of deletions overlapping one CDS")
filter(del_test, !is.na(Genome)) %>% group_by(ID, Genome) %>% count() %>% summary()
#on median average = 1, mean average =1.8

## to test the map ratio, pick out cds that will be easy to verify:
filter(del_test, !is.na(Genome)) %>% group_by(ID, Genome) %>% count() %>%
  filter(n == 1) %>% ungroup() %>% select(ID) %>% unique() %>% head()

filter(del_test, str_detect(ID, pattern = paste(c("Sobic.001G000200.1.v3.1.CDS.9",
                                            "Sobic.001G000400.5.v3.1.CDS.17",
                                            "Sobic.001G000700.2.v3.1.CDS.3",
                                            "Sobic.001G000800.1.v3.1.CDS.10",
                                            "Sobic.001G001800.1.v3.1.CDS.1",
                                            "Sobic.001G001800.1.v3.1.CDS.2"),
                                            collapse = "|"))) %>%
  select(ID, Genome, contains("Length")) %>%
  mutate(proportion_remain = 1-(Length_overlap / Length_sb)) %>% View()
#grep "Sobic.001G000200.1.v3.1.CDS.9" ZmB73_allChr.bin1_CDS.crossmap.tsv
  #map_ratio - 0.9655
#grep "Sobic.001G000400.5.v3.1.CDS.17" ZmB73_allChr.bin1_CDS.crossmap.tsv
  #map ratio = 0.9783
#grep "Sobic.001G000700.2.v3.1.CDS.3" ZmB73_allChr.bin1_CDS.crossmap.tsv
  #map ratio = 1.0000 (X)
#grep "Sobic.001G000800.1.v3.1.CDS.10" ZmB73_allChr.bin1_CDS.crossmap.tsv
  #0.9940
#grep "Sobic.001G001800.1.v3.1.CDS.1" ZmB73_allChr.bin1_CDS.crossmap.tsv
  # NA
#grep "Sobic.001G001800.1.v3.1.CDS.2" ZmB73_allChr.bin1_CDS.crossmap.tsv
  # NA


#2.
filter(del_test, !is.na(Genome)) %>%
  ggplot(aes(x=Genome, y=Length_overlap))+
  geom_boxplot(aes(color=Genome))+
  theme_bw()+xlab("")+ylab("Length of individual deletion overlaps (bp)")
filter(del_test, !is.na(Genome)) %>% select(Length_overlap) %>% summary()
#median = 8, mean = 94

filter(del_test, !is.na(Genome)) %>%
  mutate(overlap_size_class = case_when(Length_overlap <= 10 ~ "1-10bp",
                                        Length_overlap > 10 & Length_overlap <=100 ~ "11-100bp",
                                        Length_overlap > 100 & Length_overlap <=1000 ~ "101-1000bp",
                                        Length_overlap > 1000~ ">1001bp"),
         overlap_size_class = factor(overlap_size_class, levels = c("1-10bp","11-100bp","101-1000bp",">1001bp"))) %>%
  ggplot(aes(x=overlap_size_class, fill=Genome))+
  geom_bar(position = "dodge")+
  theme_bw()+xlab("Deletion overlap sizes")+ylab("Counts")
ggsave("Del_test/DelOverlapSizes_barchart.png",device="png",dpi=300)


#3. 
filter(del_test, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>%
  group_by(Genome, overlap_mult_3) %>% count() %>%
  mutate(percent = (n/69269)*100) %>%
  ggplot(aes(x=Genome, y=percent, fill=overlap_mult_3))+
  geom_bar(position="fill",stat = "identity")+
  theme_bw()+xlab("")+ylab("Proportion")+scale_fill_manual(name="Del. overlap is \nmultiple of 3bp",
                                                           values=c("grey","black"))
ggsave("Del_test/DelOverlapMult3_TFbars.png",device="png",dpi=300)
#there's a difference in the graph values and the percent values because the percent denominator counts CDS that don't have any deletions

#3.a How many CDS only have multiple of 3 deletions vs. mix vs. just non-mult. of 3 deletions?
filter(del_test, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>% 
  select(ID, Genome, overlap_mult_3) %>% unique() %>% 
  group_by(ID, Genome) %>%
  reframe(ID = ID, Genome = Genome, 
          Deletion_Types = case_when(all(isTRUE(overlap_mult_3)) ~ "Only mult 3",
                                     all(isFALSE(overlap_mult_3)) ~ "Only non-mult 3",
                                     .default = "Mixed")) %>%
  ggplot(aes(x=Deletion_Types, fill=Genome))+
  geom_bar(stat = "count", position="dodge")

#4: How can we see the congruence of map-ratio? 
B73_M1_crossmap<-read_tsv("CrossMap/NAM_cds_bed/ZmB73_allChr.bin1_CDS.crossmap.tsv",
                          col_names = c("CHROM_Zm","Start_Zm","Stop_Zm","ID","QUAL","Strand","Map_ratio"))
B73_M1_crossmap <- mutate(B73_M1_crossmap,Map_ratio = str_remove(Map_ratio, "map_ratio=") %>% as.numeric())

filter(del_test, Genome == "ZmB73") %>%
  select(ID, Genome, Length_overlap,Length_sb) %>%
  group_by(ID) %>% 
  reframe(ID=ID, Length_sb=Length_sb, total_overlap=sum(Length_overlap)) %>%
  unique() %>%
  mutate(proportion_remain = 1-(total_overlap / Length_sb)) %>%
  inner_join(y=select(B73_M1_crossmap,ID,Map_ratio), by="ID") %>%
  mutate(propRetained_diff = proportion_remain - Map_ratio) %>%
  ggplot(aes(x=propRetained_diff, y="B73"))+
  geom_jitter(width = 0)+
  geom_vline(xintercept = 0)
#negative means that the cross-map ratio is much bigger than the deletion overlap from GVCF
#mean = -0.01, median = 0 

#4.1 Are the big disagreements due to large indels in the GVCF calls?
IDs_ofInterest<-filter(del_test, Genome == "ZmB73") %>%
  select(ID, Genome, Length_overlap,Length_sb) %>%
  group_by(ID) %>% 
  reframe(ID=ID, Length_sb=Length_sb, total_overlap=sum(Length_overlap)) %>%
  unique() %>%
  mutate(proportion_remain = 1-(total_overlap / Length_sb)) %>%
  inner_join(y=select(B73_M1_crossmap,ID,Map_ratio), by="ID") %>%
  mutate(propRetained_diff = proportion_remain - Map_ratio) %>% filter(propRetained_diff < -0.8) %>% pull(ID)

filter(del_test, ID %in% IDs_ofInterest) %>% 
  mutate(del_length_size = case_when(Length_del <= 10 ~ "1-10bp lost",
                                     Length_del > 10 & Length_del <= 100 ~ "11-100bp lost",
                                     Length_del > 100 & Length_del <= 1000 ~ "101-1000bp lost",
                                     Length_del > 1000 & Length_del <= 10000 ~ "1K-10Kbp lost",
                                     Length_del > 10000 ~ ">10Kbp lost"),
         del_length_size = factor(del_length_size, levels = c("1-10bp lost","11-100bp lost","101-1000bp lost","1K-10Kbp lost",">10Kbp lost"))) %>%
  ggplot(aes(x=del_length_size, fill=Genome))+
  geom_bar(stat="count",position="dodge")

mutate(del_test, del_length_size = case_when(Length_del <= 10 ~ "1-10bp lost",
                                   Length_del > 10 & Length_del <= 100 ~ "11-100bp lost",
                                   Length_del > 100 & Length_del <= 1000 ~ "101-1000bp lost",
                                   Length_del > 1000 & Length_del <= 10000 ~ "1K-10Kbp lost",
                                   Length_del > 10000 ~ ">10Kbp lost"),
       del_length_size = factor(del_length_size, levels = c("1-10bp lost","11-100bp lost","101-1000bp lost","1K-10Kbp lost",">10Kbp lost"))) %>%
  ggplot(aes(x=del_length_size, fill=Genome))+
  geom_bar(stat="count",position="dodge")
#pretty sure this shows that these IDs that have the most disagreement are the ones that have the really big deletions around them from the GVCF

#5: Try to add in fractionation status and multiples of 3 deletion overlaps
full.fractionation.status

mult_3_fract_test<-filter(del_test, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>% 
  select(ID, Genome, overlap_mult_3) %>% unique() %>% 
  group_by(ID, Genome) %>%
  reframe(ID = ID, Genome = Genome, 
          Deletion_Types = case_when(all(isTRUE(overlap_mult_3)) ~ "Only mult 3",
                                     all(isFALSE(overlap_mult_3)) ~ "Only non-mult 3",
                                     .default = "Mixed")) %>% 
  unique()

non_reframe_test<-filter(del_test, !is.na(Genome)) %>%
  mutate(overlap_mult_3 = Length_overlap %% 3 == 0) %>% 
  select(ID, Genome, overlap_mult_3) %>% unique() %>% 
  group_by(ID, Genome) %>%
  mutate(Deletion_Types = case_when(all(isTRUE(overlap_mult_3)) ~ "Only mult 3",
                                    all(isFALSE(overlap_mult_3)) ~ "Only non-mult 3",
                                    .default = "Mixed")) %>% 
  select(-overlap_mult_3) %>%
  unique() 

mult_3_fract_test<- pivot_longer(mult_3_fract_test, cols=ends_with(".M1"), 
             names_to = "Fract_Genome", 
             values_to = "Fractionation_Call") %>%
  mutate(Fract_Genome = str_remove(Fract_Genome, ".M1")) %>%
  filter(Genome == Fract_Genome) %>% select(-Fract_Genome) 

mult_3_fract_test %>% group_by(Deletion_Types, Genome) %>% count(Fractionation_Call)
#they'll all be called fractionated because they all overlap with deletions
mult_3_fract_test %>% 
  mutate(Fractionation_Revision = case_when(Deletion_Types == "Only mult 3" ~ 0,
                                            Deletion_Types != "Only mult 3" ~ 1)) %>%
  group_by(Genome) %>% count(Fractionation_Revision)

mult_3_fract_test %>% select(-Fractionation_Call) %>% 
  mutate(gene_ID = str_split(ID, ";", simplify=TRUE)[,2] %>% str_remove("Parent=")) %>%
  group_by(gene_ID, Genome) %>% 
  reframe(gene_Deletion_Types = if(length(unique(Deletion_Types)) > 1){"Mixed"}else{unique(Deletion_Types)}) %>%
  write_tsv("Del_test/mult_3_fract_test.tsv")
  
M1_dels_by_gene<-read_tsv("M1_dels_by_gene.tsv")
filter(M1_dels_by_gene, is.na(gene_size))

M1_dels_by_gene %>% group_by(Genome, gene_size) %>%
  summarize(mean.total_overlap_del_percent = mean(total_overlap_del_percent, na.rm = T),
            sd.total_overlap_del_percent = sd(total_overlap_del_percent, na.rm=T)) %>%
  inner_join(y=count(group_by(M1_dels_by_gene, Genome, gene_size))) %>%
  mutate(se.total_overlap_del_percent = sd.total_overlap_del_percent/sqrt(n),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                          "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                          "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001"
                                          ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
         gene_size = factor(gene_size, levels = c("<= 1Kbp","1-2Kbp","2-3Kbp","3-4Kbp","4-5Kbp","> 5Kbp")))%>%
  ggplot(aes(x=gene_size, y=mean.total_overlap_del_percent))+
  geom_bar(aes(fill=Genome), stat = "identity",position = "dodge")+
  geom_errorbar(aes(ymin=mean.total_overlap_del_percent-se.total_overlap_del_percent,
                    ymax=mean.total_overlap_del_percent+se.total_overlap_del_percent,
                    color=Genome), position = "dodge")+
  scale_fill_manual(values=genome_colors)+ scale_color_manual(values = rep("#000000",36))+
  guides(color="none", 
         fill = guide_legend(ncol=10))+ 
  theme_bw()+ theme(legend.position = "bottom")+
  xlab("Gene length binned")+ylab("Mean Percent Gene \nCoding Sequence Deleted")
ggsave("../Fractionation_Plots/M1_GeneLengthBins_by_PercentGeneDelOverlap.png",
       device = "png", dpi=300, width = 4.5, height = 4.5, units = "in")

M2_dels_by_gene<-read_tsv("M2_dels_by_gene.tsv")

M2_dels_by_gene %>% group_by(Genome, gene_size) %>%
  summarize(mean.total_overlap_del_percent = mean(total_overlap_del_percent, na.rm = T),
            sd.total_overlap_del_percent = sd(total_overlap_del_percent, na.rm=T)) %>%
  inner_join(y=count(group_by(M1_dels_by_gene, Genome, gene_size))) %>%
  mutate(se.total_overlap_del_percent = sd.total_overlap_del_percent/sqrt(n),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                          "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                          "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001"
                                          ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
         gene_size = factor(gene_size, levels = c("<= 1Kbp","1-2Kbp","2-3Kbp","3-4Kbp","4-5Kbp","> 5Kbp")))%>%
  ggplot(aes(x=gene_size, y=mean.total_overlap_del_percent))+
  geom_bar(aes(fill=Genome), stat = "identity",position = "dodge")+
  geom_errorbar(aes(ymin=mean.total_overlap_del_percent-se.total_overlap_del_percent,
                    ymax=mean.total_overlap_del_percent+se.total_overlap_del_percent,
                    color=Genome), position = "dodge")+
  scale_fill_manual(values=genome_colors)+ scale_color_manual(values = rep("#000000",36))+
  guides(color="none", 
         fill = guide_legend(ncol=10))+ 
  theme_bw()+ theme(legend.position = "bottom")+
  xlab("Gene length binned")+ylab("Mean Percent Gene \nCoding Sequence Deleted")
ggsave("../Fractionation_Plots/M2_GeneLengthBins_by_PercentGeneDelOverlap.png",
       device = "png", dpi=300, width = 4.5, height = 4.5, units = "in")

taxon_colors<-c("maize, temperate"= "#92ddb0",
                "maize, mixed"= "#7ccbb2",
                "maize, tropical"= "#4da89d",
                "parviglumis"= "#018ce7",
                "mexicana"= "#036db2",
                "huehuetenangensis"= "#2B5A78",
                "diploperennis"= "#433475",
                "tripsacum"= "#F45B69")

dels_by_gene<-add_column(M1_dels_by_gene, Subgenome = "M1") %>% add_row(M2_dels_by_gene, Subgenome = "M2")
dels_by_gene %>% group_by(Genome, Subgenome) %>%
  summarize(mean.total_overlap_del_percent = mean(total_overlap_del_percent, na.rm = T),
            sd.total_overlap_del_percent = sd(total_overlap_del_percent, na.rm=T)) %>%
  inner_join(y=count(group_by(dels_by_gene, Genome, Subgenome))) %>%
  mutate(se.total_overlap_del_percent = sd.total_overlap_del_percent/sqrt(n),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                          "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                          "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001"
                                          ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
         Taxon = case_when(Genome %in% c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                         "ZmHP301") ~ "maize, temperate",
                           Genome %in% c("ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ "maize, mixed",
                           Genome %in% c("ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                         "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                         "ZmNC350" ,"ZmNC358" ,"ZmTzi8") ~ "maize, tropical",
                           str_starts(Genome, "Zv") ~ "parviglumis",
                           str_starts(Genome, "Zx") ~ "mexicana",
                           str_starts(Genome, "Zh") ~ "huehuetenangensis",
                           str_starts(Genome, "Zd") ~ "diploperennis",
                           str_starts(Genome, "Td") ~ "tripsacum") %>% factor(levels = c("maize, temperate","maize, mixed","maize, tropical", "parviglumis","mexicana","huehuetenangensis","diploperennis","tripsacum")),
         )%>%
  ggplot(aes(x=Subgenome, y=mean.total_overlap_del_percent))+
  geom_bar(aes(fill=Genome), stat = "identity",position = "dodge")+
  geom_errorbar(aes(ymin=mean.total_overlap_del_percent-se.total_overlap_del_percent,
                    ymax=mean.total_overlap_del_percent+se.total_overlap_del_percent,
                    color=Genome), position = "dodge")+
  scale_fill_manual(values=genome_colors)+ scale_color_manual(values = rep("#000000",36))+
  guides(color="none", 
         fill = "none")+ 
  theme_bw()+ 
  xlab("Subgenome")+ylab("Mean Percent Gene \nCoding Sequence Deleted")
ggsave("../Fractionation_Plots/BothSubgenomes_GeneLengthBins_by_PercentGeneDelOverlap.png",
       device = "png", dpi=300, width = 4.5, height = 4.5, units = "in")

all_dels<-read_tsv("intermediate-data-files/M1_mult_3_deletions.tsv") %>% 
  add_column(Subgenome = "M1") %>%
  add_row(read_tsv("intermediate-data-files/M2_mult_3_deletions.tsv"), Subgenome = "M2")

filter(all_dels, !is.na(Genome)) %>%
  mutate(frame_shift = case_when(Deletion_Types %in% c("Mixed", "Only non-mult 3") ~ T,
                                 Deletion_Types == "Only mult 3" ~ F),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                          "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                          "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001"
                                          ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS"))
  ) %>% 
  select(ID, Genome, Subgenome, frame_shift) %>% unique() %>% 
  ggplot(aes(x=frame_shift, fill=Genome))+
  geom_bar(stat = "count", position="dodge")+
  guides(fill = "none")+scale_fill_manual(values = genome_colors)+
  facet_wrap(vars(Subgenome))+theme_bw()+
  xlab("Contain Frame-Shifting Deletion")+
  ylab("Counts of Exons with Deletions")
ggsave("../Fractionation_Plots/BothSubgenomes_Mult3DelTypes_barchart.png",
       device = "png", dpi=300, width = 4.5, height = 4.5, units = "in")

all_dels<-all_dels %>% 
  mutate(frame_shift = case_when(Deletion_Types %in% c("Mixed", "Only non-mult 3") ~ T,
                               Deletion_Types == "Only mult 3" ~ F),
       Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                        "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                        "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,
                                        "ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,
                                        "ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                        "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001"
                                        ,"ZdGigi_4to1","ZdMomo_4to1","TdFL","TdKS")),
       )

delType_byStatus<-full_join(long.full.fractionation.status, all_dels, by=c("ID", "Genome"))
delType_byStatus<-filter(delType_byStatus, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697","ZnPI615697_4to1"))

delType_byStatus<- delType_byStatus  %>%
  mutate(Subgenome = case_when(is.na(Subgenome) & Status == "M1_Retained:M2_NA" ~ "M1",
                               is.na(Subgenome) & Status == "M1_NA:M2_Retained" ~ "M2",
                               .default = Subgenome))
delType_byStatus<-delType_byStatus %>% 
  mutate(Subgenome = case_when(is.na(Subgenome) & Status == "Both_Retained" ~ "M1",
                               .default = Subgenome))
#got to add in the M2 row for "Both_Retained"
delType_byStatus<-delType_byStatus %>% filter(Status == "Both_Retained") %>% unique() %>%
  mutate(Subgenome = "M2") %>%
  add_row(delType_byStatus) %>% unique()

delType_byStatus<-mutate(delType_byStatus, 
      Overall_Status = case_when(frame_shift==TRUE ~ "Frame-shift Fractionation",
                                  frame_shift == FALSE ~ "Non-frameshift Fractionation",
                                  is.na(frame_shift) & Status == "Both_NA" ~ NA,
                                  is.na(frame_shift) & Status %in% c("Both_Retained", "M2_Retained","M1_Retained","M1_Retained:M2_NA","M1_NA:M2_Retained") ~ "Retained"),
       Overall_Status = factor(Overall_Status, levels = c("Retained","Non-frameshift Fractionation","Frame-shift Fractionation")))

delType_byStatus %>% filter(Status != "Both_NA") %>%
  group_by(Genome, Subgenome, Overall_Status) %>% 
  count() %>%
  mutate(Percentage = (n/69269)*100) %>%
  ggplot(aes(x=Overall_Status, y=Percentage))+
  geom_violin(aes(fill = Subgenome))+
  theme_bw()+xlab("")+ylab("Percentage of Exons")+
  scale_fill_manual(values = subgenome_colors)
ggsave("../Fractionation_Plots/delType_byStatus.violin.bySubgenome.png",
       device = "png", dpi=300)

delType_byStatus %>% filter(Status != "Both_NA") %>%
  group_by(Genome, Subgenome, Overall_Status) %>% 
  count() %>%
  mutate(Percentage = (n/69269)*100) %>%
  ggplot(aes(x=Overall_Status, y=Percentage))+
  geom_jitter(aes(color = Genome), height = 0, width = 0.25, alpha =0.75)+
  geom_violin(alpha = 0)+
  theme_bw()+xlab("")+ylab("Percentage of Exons")+
  scale_color_manual(values = genome_colors)+
  facet_grid(cols = vars(Subgenome))+
  guides(color = "none")+
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))+
  scale_x_discrete(labels = c("Retained","Non-frameshift\nFractionation","Frame-shift\nFractionation"))
ggsave("../Fractionation_Plots/delType_byStatus.violin.byGenome.png",
       device = "png", dpi=300, height = 5, width = 7, units = "in")

delType_byStatus %>% filter(Status != "Both_NA") %>%
  group_by(Genome, Subgenome, Overall_Status) %>% 
  count() %>%
  mutate(Percentage = (n/69269)*100) %>%
  group_by(Subgenome, Overall_Status) %>%
  summarize(mean.percent = mean(Percentage,na.rm=T))


#Dn/Ds Convergence
convergence.NoZnZd<-read_tsv("convergence.NoZnZd.tsv")
convergence.dNdS<- convergence.NoZnZd %>%
  mutate(Gene_ID = str_remove_all(Gene_ID, ".1.v3.1")) %>%
  filter(Gene_ID %in% dn_ds_data$GeneID_Sb313)

#Grab just the B73 lines
convergence.sharing<-convergence.dNdS %>% 
  filter(Genome == 'ZmB73') 
#grab the other genomes and add to create a wider dataframe
convergence.sharing<-filter(convergence.dNdS, Genome != "ZmB73") %>%
  left_join(convergence.sharing, by=c("Gene_ID","M"), 
            suffix = c(".other",".B73"))

#Make the pattern matching :)
convergence.sharing<-convergence.sharing %>%
  mutate(Sharing = case_when(Loss_Pattern.B73 == Loss_Pattern.other ~ "Complete_Sharing"))

for(i in 1:nrow(convergence.sharing)){
#for(i in 1:100){ #to test loops
  if(is.na(convergence.sharing$Sharing[i])){
    if(is.na(convergence.sharing$Loss_Pattern.B73[i]) | is.na(convergence.sharing$Loss_Pattern.other[i])){
      convergence.sharing$Sharing[i]<-"Retention_Difference"
    }else{
    loss_pattern_B73<-str_split(convergence.sharing$Loss_Pattern.B73[i], ":", simplify=T)
    loss_pattern_other<-str_split(convergence.sharing$Loss_Pattern.other[i],":", simplify = T)
    
    if(length(loss_pattern_other) > length(loss_pattern_B73)){
      long_pattern<-loss_pattern_other
      short_pattern<-loss_pattern_B73
    }else{
      long_pattern<-loss_pattern_B73
      short_pattern<-loss_pattern_other
    }
    #count the number of exons that are detected from shortest pattern in longest pattern
    #count the number of exons that are detected from shortest pattern in longest pattern
    counter<-0
    for(e in 1:length(short_pattern)){
      if(any(str_detect(long_pattern, short_pattern[1,e]))){
        counter<-counter+1
      }
    }
    if(counter > 0){
      convergence.sharing$Sharing[i]<-"Some_Sharing"
    }else{
      convergence.sharing$Sharing[i]<-"Completely_Different"
    }
  }
  }
}

convergence.sharing %>% group_by(Sharing) %>% count()

#add in the dN/dS values

dnds<-Yin2022MBE_data %>% 
  select("ω_M1","ω_M2", "GeneID_Sb313") %>%
  pivot_longer(cols=c("ω_M1","ω_M2"),names_to = "Subgenome",values_to = "dNdS") %>%
  mutate(Subgenome = str_replace(Subgenome, "ω_","")) %>%
  right_join(convergence.sharing, by = c("GeneID_Sb313"="Gene_ID", "Subgenome"="M"))

#NEXT STEP: Plot dN/dS values by Sharing

dnds<-mutate(dnds, Sharing = case_when(Sharing == "Complete_Sharing" ~ "Frac.: Complete Share",
                                 Sharing == "Completely_Different" ~ "Frac.: No Share",
                                 Sharing == "Some_Sharing" ~ "Frac.: Some Share",
                                 Sharing == "Retention_Difference" ~ "Ret.: No Share") %>% factor(levels = c("Frac.: Complete Share","Frac.: Some Share","Frac.: No Share","Ret.: No Share")))
dnds<-dnds %>% 
  mutate(Genomes_Compared = case_when(
    str_detect(Genome.other, "Zm") ~ "Zea_accessions",
    str_detect(Genome.other, "Zv") ~ "Zea_subspecies",
    str_detect(Genome.other, "Zx") ~ "Zea_subspecies",
    str_detect(Genome.other, "Zh") ~ "Zea_subspecies",
    str_detect(Genome.other, "Zd") ~ "Zea_species",
    str_detect(Genome.other, "Td") ~ "Trip_v_Zea",
    ) %>% factor(levels = c("Zea_accessions","Zea_subspecies","Zea_species","Trip_v_Zea"))
  )

#1. What is overall ZmB73 dN/dS of this subset? 
dn_ds_data %>% summarize(mean.M1 = mean(ω_M1),
                         mean.M2 = mean(ω_M2),
                         med.M1 = median(ω_M1),
                         med.M2 = median(ω_M2))

t.test(x=dn_ds_data$ω_M1, y=dn_ds_data$ω_M2, paired = T)
#df = 3194, p-value = 0.04794, mean difference = -0.00933734

dn_ds_data %>% select(ω_M1,ω_M2) %>%
  pivot_longer(cols = everything(), names_to = "Subgenome", values_to = "dNdS") %>%
  mutate(Subgenome = str_replace(Subgenome, "ω_","")) %>%
  ggplot(aes(x=Subgenome, y=dNdS))+
  geom_jitter(height = 0, alpha =0.25)+
  geom_boxplot(aes(fill = Subgenome), alpha = 0.95)+
  theme_bw()+guides(fill = "none")+scale_fill_manual(values = subgenome_colors)+
  xlab("")+ylab("Zm B73 dN/dS")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))
ggsave("../Fractionation_Plots/2026-03-16-B73dNdS_Overall.pdf",
       dpi=300, device = "pdf")

#2. How / does dN/dS vary by sharing of fractionation with B73?

M1.FCS<-filter(dnds, Sharing == "Frac.: Complete Share", Subgenome == "M1") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M1.FCS) #no differences

M2.FCS<-filter(dnds, Sharing == "Frac.: Complete Share", Subgenome == "M2") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M2.FCS) #no differences

M1.FSS<-filter(dnds, Sharing == "Frac.: Some Share", Subgenome == "M1") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M1.FSS) #no differences

M2.FSS<-filter(dnds, Sharing == "Frac.: Some Share", Subgenome == "M2") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M2.FSS) #p = 7.67e-06
TukeyHSD(M2.FSS) #Trip_v_Zea vs. the other three categories

M1.FNS<-filter(dnds, Sharing == "Frac.: No Share", Subgenome == "M1") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M1.FNS) #no differences

M2.FNS<-filter(dnds, Sharing == "Frac.: No Share", Subgenome == "M2") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M2.FNS) #p=0.0133
TukeyHSD(M2.FNS) #Zea subsp vs. Zea accessions

M1.RNS<-filter(dnds, Sharing == "Ret.: No Share", Subgenome == "M1") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M1.RNS) #no differences

M2.RNS<-filter(dnds, Sharing == "Ret.: No Share", Subgenome == "M2") %>%
  aov(data = ., formula = dNdS ~ Genomes_Compared)
summary(M2.RNS) #no differences

ggplot(dnds, aes(x= Sharing, y=dNdS))+
  geom_boxplot(aes(fill = Genomes_Compared), alpha = 0.75)+
  theme_bw()+ylab("Zm B73 dN/dS")+
  xlab("")+facet_grid(vars(Subgenome))+
  scale_fill_manual(values = c("Trip_v_Zea" = "#F45B69","Trip_accession" = "#f0a3aa","Zea_species" = "#412151","Zea_subspecies" = "#247590","Zea_accessions" = "#4da89d"),
                    name = "Genomes Compared")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16))
ggsave('../Fractionation_Plots/2026-03-16-dNdS_bySharing_byGenomesCompared.pdf',
       dpi = 300, device = "pdf", width = 8)

## No differences when we split out how diverged genomes are in the genomes compared

dnds.aov<-aov(dNdS ~ Subgenome+Sharing, data = dnds)
summary(dnds.aov)
TukeyHSD(dnds.aov, which = "Sharing")
#FCS = B, FSS = A, FNS=B, RNS = C

ggplot(dnds, aes(x= Sharing, y=dNdS))+
  geom_boxplot(aes(fill = Subgenome))+
  theme_bw()+ylab("Zm B73 dN/dS")+
  xlab("")+scale_fill_manual(values = subgenome_colors)+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size=16))
ggsave("../Fractionation_Plots/2026-03-16-dNdS_bySharing_bySubgenome.pdf",
       device = "pdf",dpi=300,width = 8)

#3. Are patterns in 2 due to differences in sample size?

dnds.samplesizes<-dnds %>% group_by(Sharing, Subgenome) %>% count()

observed_dnds.avg<-dnds %>% group_by(Sharing, Subgenome) %>% summarize(avg.dnds = mean(dNdS, na.rm =T))

dnds_permutation<-function(){
  permutation<-tibble(dNdS = sample(dnds$dNdS),
                      Sharing = c(rep(dnds.samplesizes$Sharing[1], dnds.samplesizes$n[1]),
                                  rep(dnds.samplesizes$Sharing[2], dnds.samplesizes$n[2]),
                                  rep(dnds.samplesizes$Sharing[3], dnds.samplesizes$n[3]),
                                  rep(dnds.samplesizes$Sharing[4], dnds.samplesizes$n[4]),
                                  rep(dnds.samplesizes$Sharing[5], dnds.samplesizes$n[5]),
                                  rep(dnds.samplesizes$Sharing[6], dnds.samplesizes$n[6]),
                                  rep(dnds.samplesizes$Sharing[7], dnds.samplesizes$n[7]),
                                  rep(dnds.samplesizes$Sharing[8], dnds.samplesizes$n[8])),
                      Subgenome = c(rep(dnds.samplesizes$Subgenome[1], dnds.samplesizes$n[1]),
                                    rep(dnds.samplesizes$Subgenome[2], dnds.samplesizes$n[2]),
                                    rep(dnds.samplesizes$Subgenome[3], dnds.samplesizes$n[3]),
                                    rep(dnds.samplesizes$Subgenome[4], dnds.samplesizes$n[4]),
                                    rep(dnds.samplesizes$Subgenome[5], dnds.samplesizes$n[5]),
                                    rep(dnds.samplesizes$Subgenome[6], dnds.samplesizes$n[6]),
                                    rep(dnds.samplesizes$Subgenome[7], dnds.samplesizes$n[7]),
                                    rep(dnds.samplesizes$Subgenome[8], dnds.samplesizes$n[8])))
  averages<-permutation %>% group_by(Sharing, Subgenome) %>% summarize(avg.dnds = mean(dNdS, na.rm =T))
  return(averages)
}

dnds_permutation.test<-tibble(dnds_permutation(), Set = 1)

for(i in 2:1000){
  dnds_permutation.test<-add_row(dnds_permutation.test, 
                                 dnds_permutation(), 
                                 Set = i)
}

dnds_permutation.test<-left_join(dnds_permutation.test, 
          observed_dnds.avg, 
          by = c("Sharing","Subgenome"), 
          suffix = c(".perm",".obs"))

ggplot(dnds_permutation.test, aes(x=Sharing))+
  geom_boxplot(aes(fill = Subgenome, y=avg.dnds.perm))+
  geom_point(aes(y=avg.dnds.obs, color = Subgenome),
             position = position_dodge(width = 0.75))+
  theme_bw()+xlab("")+ylab("Permutation Mean Zm B73 dN/dS")+
  scale_fill_manual(values = subgenome_colors)+
  scale_color_manual(values = c("M1"="#0055ff", "M2"="#db021f"))+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size=16))
ggsave("../Fractionation_Plots/2026-03-16-dNdS-permutationResults.pdf",
       dpi=300, device = "pdf", width = 8)

dnds_permutation.test %>% group_by(Sharing, Subgenome) %>%
  reframe(mean.perm =mean(avg.dnds.perm), sd.perm = sd(avg.dnds.perm)) %>%
  add_column(avg.dnds.obs = c(0.254,0.269,0.347,0.391,0.234,0.306,0.205,0.228)) %>%
  mutate(lowerCI = mean.perm - 3*sd.perm,
         upperCI = mean.perm + 3*sd.perm,
         ObsInCI = case_when(avg.dnds.obs <= upperCI & avg.dnds.obs >= lowerCI ~ "YES",
                             .default = "NO"))


## Differential GO: gene length and exon number question
differential_retention.NoZnZd<-mutate(differential_retention.NoZnZd,
                                      Status_number = case_when(Pattern %in% c("BothDeleted","M2Retained","M1Retained","BothRetained") ~ 1, 
                                                                Pattern %in% c("BothRetained:M2Retained","BothDeleted:M2Retained","BothDeleted:M1Retained","BothRetained:M1Retained","M1Retained:M2Retained","BothDeleted:BothRetained") ~ 2,
                                                                Pattern %in% c("BothRetained:M2Retained:BothDeleted","BothRetained:M1Retained:BothDeleted","M1Retained:M2Retained:BothDeleted","BothRetained:M1Retained:M2Retained") ~ 3,
                                                                Pattern %in% c("BothRetained:M1Retained:M2Retained:BothDeleted") ~ 4,
                                                                .default = NA)
                                      )

temp_ref.gene<-mutate(ref_Sb313.cds,
                      Gene_ID = str_split(ID, ";", simplify = T)[,2] %>% str_remove("Parent="))
temp_ref.gene<-temp_ref.gene %>% 
  group_by(Gene_ID) %>%
  reframe(gene_start = min(Start), 
          gene_stop = max(Stop),
          Chr = Chr,
          exonCt = n()) %>% unique()

DFGO_tests<-differential_retention.NoZnZd %>% 
  select(Gene_ID, Status_number) %>% 
  unique() %>% na.omit() %>%
  left_join(temp_ref.gene)

DFGO_tests<-mutate(DFGO_tests, 
                   Gene_Length = gene_stop - gene_start,
                   Status_number = factor(Status_number, levels = c(1,2,3,4)))

ggplot(DFGO_tests, aes(x=Status_number, y=Gene_Length))+
  geom_boxplot()
aov(Gene_Length ~ Status_number, data = DFGO_tests) %>% summary()
TukeyHSD(aov(Gene_Length ~ Status_number, data = DFGO_tests))
#Single status genes are statistically shorter than multiple status genes by 446-630bps

ggplot(DFGO_tests, aes(x=Status_number, y=exonCt))+
  geom_boxplot()
aov(exonCt ~ Status_number, data = DFGO_tests) %>% summary()
TukeyHSD(aov(exonCt ~ Status_number, data = DFGO_tests))
#Single status genes have statistically 1-1.3 fewer exons than all other statuses

#https://pmc.ncbi.nlm.nih.gov/articles/PMC137605/1

# table counts by taxon
long.full.fractionation.status %>% 
  filter(!Genome %in% c("ZdGigi","ZdMomo")) %>%
  mutate(Taxon = case_when(str_starts(Genome, "Zm") ~ "maize",
                           str_starts(Genome, "Zv") ~ "parviglumis",
                    str_starts(Genome, "Zx") ~ "mexicana",
                    str_starts(Genome, "Zh") ~ "huehuetenangensis",
                    str_starts(Genome, "Zd") ~ "diploperennis",
                    str_starts(Genome, "Zn") ~ "nicaraguensis",
                    str_starts(Genome, "Td") ~ "tripsacum")) %>%
  group_by(Taxon, Status) %>%
  count() %>%
  mutate(avg.count = case_when(Taxon == "diploperennis" ~ n/2,
                               Taxon == "tripsacum" ~ n/2,
                               Taxon == "nicaraguensis" ~ n,
                               Taxon == "huehuetenangensis" ~ n,
                               Taxon == "mexicana" ~ n/2,
                               Taxon == "parviglumis" ~ n/2,
                               Taxon == "maize" ~ n/26)) %>%
  View()
