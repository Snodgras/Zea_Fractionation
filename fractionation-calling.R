#To look at deletion calls as fractionation events using the full set of genomes and Chromosomes
install.packages("multcomp")
library(tidyverse)

####MAKING THE REFERENCE CDS COORDINATES####
#The specific gene models to include as reference Sb313
Sb313_ref_models<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt", col_names = c("Sb313_gene_model","Sb_exon_count","Av_gene_model","Av_exon_count","Proportion_AvExonCnt_By_SbExonCnt"))
#reduce isoforms
##make a column that's just the gene ID part of the model name
Sb313_ref_models<-mutate(Sb313_ref_models, Sb313_gene_ID = str_remove(Sb313_gene_model,pattern ="\\.[0-9].v3.1"))
Sb313_ref_models$Sb313_gene_ID<-str_remove(Sb313_ref_models$Sb313_gene_ID, pattern = "\\.[0-9][0-9].v3.1")
uniq_ref_Sb313_gene_ID<-unique(Sb313_ref_models$Sb313_gene_ID)

##remove isoforms such that keep the one with most exons and if a tie, pick one at random
clean_Sb313_ref_models<-tibble(Sb313_gene_model = NA, 
                               Sb_exon_count=NA,
                               Av_gene_model=NA,
                               Av_exon_count=NA,
                               Proportion_AvExonCnt_By_SbExonCnt=NA,
                               Sb313_gene_ID=NA)
for(i in 1:length(uniq_ref_Sb313_gene_ID)){
  t<-filter(Sb313_ref_models, Sb313_gene_ID == uniq_ref_Sb313_gene_ID[i])
  if(nrow(t) > 1){ #if there's isoforms
    u<-grep(x=t$Sb_exon_count, pattern=max(t$Sb_exon_count)) #find which row has the most exons
    if(length(u) > 1){ #if multiple isoforms have the max exon count
      v<-sample(u,1) #randomly pick one
      clean_Sb313_ref_models<-add_row(clean_Sb313_ref_models,Sb313_gene_model = pull(t[v,1]), 
                                      Sb_exon_count=pull(t[v,2]),
                                      Av_gene_model=pull(t[v,3]),
                                      Av_exon_count=pull(t[v,4]),
                                      Proportion_AvExonCnt_By_SbExonCnt=pull(t[v,5]),
                                      Sb313_gene_ID=pull(t[v,6])) #add that one to the clean models df
    }else{ #if just 1 isoform has the most exons, add that one to the clean models df
      clean_Sb313_ref_models<-add_row(clean_Sb313_ref_models,Sb313_gene_model = pull(t[u,1]), 
                                      Sb_exon_count=pull(t[u,2]),
                                      Av_gene_model=pull(t[u,3]),
                                      Av_exon_count=pull(t[u,4]),
                                      Proportion_AvExonCnt_By_SbExonCnt=pull(t[u,5]),
                                      Sb313_gene_ID=pull(t[u,6]))
    }
  }else{#if there are no isoforms
    clean_Sb313_ref_models<-add_row(clean_Sb313_ref_models,Sb313_gene_model = pull(t[,1]), 
                                    Sb_exon_count=pull(t[,2]),
                                    Av_gene_model=pull(t[,3]),
                                    Av_exon_count=pull(t[,4]),
                                    Proportion_AvExonCnt_By_SbExonCnt=pull(t[,5]),
                                    Sb313_gene_ID=pull(t[,6]))
  }
}
clean_Sb313_ref_models<-clean_Sb313_ref_models[-1,]

#subset the CDS coordinates with the ref models
Sb313.cds<-mutate(Sb313.cds, CDS_ID = str_split(ID,pattern = ";", simplify = T) %>% .[,1],
                  Gene_ID = str_split(ID,pattern = ";", simplify = T) %>% .[,2])
Sb313.cds$Gene_ID<-str_remove(Sb313.cds$Gene_ID,pattern = "Parent=")
Sb313.cds$CDS_ID<-str_remove(Sb313.cds$CDS_ID,pattern = "ID=")

ref_Sb313.cds<-tibble(CHROM=NA, Start=NA, End=NA, ID=NA, Strand=NA, CDS_ID=NA,Gene_ID=NA)
for(i in 1:nrow(clean_Sb313_ref_models)){
  t<-filter(Sb313.cds, Gene_ID == clean_Sb313_ref_models$Sb313_gene_model[i])
  ref_Sb313.cds<-add_row(ref_Sb313.cds, CHROM=pull(t[,1]), Start=pull(t[,2]), 
                         End=pull(t[,3]), ID=pull(t[,4]), Strand=pull(t[,5]), 
                         CDS_ID=pull(t[,6]),Gene_ID=pull(t[,7]))
}
ref_Sb313.cds<-ref_Sb313.cds[-1,]

ref_Sb313.cds<-ref_Sb313.cds %>% mutate(CDS_Length = End - Start) 
ref_Sb313.cds %>% mutate(score = ".") %>% 
  dplyr::select(c(CHROM, Start, End, ID, score, Strand)) %>%
  write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed", col_names = FALSE)

####BRINGING IN THE DELETED/RETAINED/NOALIGN CALLED####
tripsacinae_genome_IDs<-c("TdFL", "TdKS", "ZdGigi_4to1", "ZdMomo_4to1", "ZnPI615697_4to1", "ZdGigi", "ZdMomo", "ZhRIMHU001", "ZmB73", "ZmB97", "ZmCML103", "ZmCML228",
                          "ZmCML247", "ZmCML277", "ZmCML322", "ZmCML333", "ZmCML52", "ZmCML69", "ZmHP301", "ZmIL14H", "ZmKi11",
                          "ZmKi3", "ZmKy21", "ZmM162W", "ZmM37W", "ZmMo18W", "ZmMS71", "ZmNC350", "ZmNC358", "ZmOh43", "ZmOh7b",
                          "ZmP39", "ZmTx303", "ZmTzi8", "ZnPI615697", "ZvTIL01", "ZvTIL11", "ZxTIL18", "ZxTIL25")
#genome.M.ID.del

for(i in 1:length(tripsacinae_genome_IDs)){
  assign(paste0(tripsacinae_genome_IDs[i],".M1",".ID.del"), 
         read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/",
                         tripsacinae_genome_IDs[i],".bin1.allchr.refExons.deleted")))
  assign(paste0(tripsacinae_genome_IDs[i],".M2",".ID.del"), 
         read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/",
                         tripsacinae_genome_IDs[i],".bin2.allchr.refExons.deleted")))
}

#genome.M.ID.ret
for(i in 1:length(tripsacinae_genome_IDs)){
  assign(paste0(tripsacinae_genome_IDs[i],".M1",".ID.ret"), 
         read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/",
                         tripsacinae_genome_IDs[i],".bin1.allchr.refExons.retained")))
  assign(paste0(tripsacinae_genome_IDs[i],".M2",".ID.ret"), 
         read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/",
                         tripsacinae_genome_IDs[i],".bin2.allchr.refExons.retained")))
}

#genome.M.ID.NA
for(i in 1:length(tripsacinae_genome_IDs)){
  assign(paste0(tripsacinae_genome_IDs[i],".M1",".ID.NA"), 
         read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/",
                         tripsacinae_genome_IDs[i],".bin1.allchr.refExons.noalign")))
  assign(paste0(tripsacinae_genome_IDs[i],".M2",".ID.NA"), 
         read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/",
                         tripsacinae_genome_IDs[i],".bin2.allchr.refExons.noalign")))
}

####JOIN CALLS ACROSS GENOMES####
full.fractionation.status<-tibble(ID = filter(ref_Sb313.cds, !str_detect(CHROM,"super")) %>% select(ID) %>% pull())
full.fractionation.status<-full.fractionation.status %>% 
  mutate(TdFL.M1 = case_when(ID %in% pull(TdFL.M1.ID.del[,1]) ~ 1, ID %in% pull(TdFL.M1.ID.ret[,1]) ~ 0, ID %in% pull(TdFL.M1.ID.NA[,1]) ~ NA),
         TdFL.M2 = case_when(ID %in% pull(TdFL.M2.ID.del[,1]) ~ 1, ID %in% pull(TdFL.M2.ID.ret[,1]) ~ 0, ID %in% pull(TdFL.M2.ID.NA[,1]) ~ NA),
         TdKS.M1 = case_when(ID %in% pull(TdKS.M1.ID.del[,1]) ~ 1, ID %in% pull(TdKS.M1.ID.ret[,1]) ~ 0, ID %in% pull(TdKS.M1.ID.NA[,1]) ~ NA),
         TdKS.M2 = case_when(ID %in% pull(TdKS.M2.ID.del[,1]) ~ 1, ID %in% pull(TdKS.M2.ID.ret[,1]) ~ 0, ID %in% pull(TdKS.M2.ID.NA[,1]) ~ NA),
         ZdGigi_4to1.M1 = case_when(ID %in% pull(ZdGigi_4to1.M1.ID.del[,1]) ~ 1, ID %in% pull(ZdGigi_4to1.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZdGigi_4to1.M1.ID.NA[,1]) ~ NA),
         ZdGigi_4to1.M2 = case_when(ID %in% pull(ZdGigi_4to1.M2.ID.del[,1]) ~ 1, ID %in% pull(ZdGigi_4to1.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZdGigi_4to1.M2.ID.NA[,1]) ~ NA),
         ZdMomo_4to1.M1 = case_when(ID %in% pull(ZdMomo_4to1.M1.ID.del[,1]) ~ 1, ID %in% pull(ZdMomo_4to1.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZdMomo_4to1.M1.ID.NA[,1]) ~ NA),
         ZdMomo_4to1.M2 = case_when(ID %in% pull(ZdMomo_4to1.M2.ID.del[,1]) ~ 1, ID %in% pull(ZdMomo_4to1.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZdMomo_4to1.M2.ID.NA[,1]) ~ NA),
         ZnPI615697_4to1.M1 = case_when(ID %in% pull(ZnPI615697_4to1.M1.ID.del[,1]) ~ 1, ID %in% pull(ZnPI615697_4to1.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZnPI615697_4to1.M1.ID.NA[,1]) ~ NA),
         ZnPI615697_4to1.M2 = case_when(ID %in% pull(ZnPI615697_4to1.M2.ID.del[,1]) ~ 1, ID %in% pull(ZnPI615697_4to1.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZnPI615697_4to1.M2.ID.NA[,1]) ~ NA),
         ZdGigi.M1 = case_when(ID %in% pull(ZdGigi.M1.ID.del[,1]) ~ 1, ID %in% pull(ZdGigi.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZdGigi.M1.ID.NA[,1]) ~ NA),
         ZdGigi.M2 = case_when(ID %in% pull(ZdGigi.M2.ID.del[,1]) ~ 1, ID %in% pull(ZdGigi.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZdGigi.M2.ID.NA[,1]) ~ NA),
         ZdMomo.M1 = case_when(ID %in% pull(ZdMomo.M1.ID.del[,1]) ~ 1, ID %in% pull(ZdMomo.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZdMomo.M1.ID.NA[,1]) ~ NA),
         ZdMomo.M2 = case_when(ID %in% pull(ZdMomo.M2.ID.del[,1]) ~ 1, ID %in% pull(ZdMomo.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZdMomo.M2.ID.NA[,1]) ~ NA),
         ZhRIMHU001.M1 = case_when(ID %in% pull(ZhRIMHU001.M1.ID.del[,1]) ~ 1, ID %in% pull(ZhRIMHU001.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZhRIMHU001.M1.ID.NA[,1]) ~ NA),
         ZhRIMHU001.M2 = case_when(ID %in% pull(ZhRIMHU001.M2.ID.del[,1]) ~ 1, ID %in% pull(ZhRIMHU001.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZhRIMHU001.M2.ID.NA[,1]) ~ NA),
         ZmB73.M1 = case_when(ID %in% pull(ZmB73.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmB73.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmB73.M1.ID.NA[,1]) ~ NA),
         ZmB73.M2 = case_when(ID %in% pull(ZmB73.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmB73.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmB73.M2.ID.NA[,1]) ~ NA),
         ZmB97.M1 = case_when(ID %in% pull(ZmB97.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmB97.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmB97.M1.ID.NA[,1]) ~ NA),
         ZmB97.M2 = case_when(ID %in% pull(ZmB97.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmB97.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmB97.M2.ID.NA[,1]) ~ NA),
         ZmCML103.M1 = case_when(ID %in% pull(ZmCML103.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML103.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML103.M1.ID.NA[,1]) ~ NA),
         ZmCML103.M2 = case_when(ID %in% pull(ZmCML103.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML103.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML103.M2.ID.NA[,1]) ~ NA),
         ZmCML228.M1 = case_when(ID %in% pull(ZmCML228.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML228.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML228.M1.ID.NA[,1]) ~ NA),
         ZmCML228.M2 = case_when(ID %in% pull(ZmCML228.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML228.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML228.M2.ID.NA[,1]) ~ NA),
         ZmCML247.M1 = case_when(ID %in% pull(ZmCML247.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML247.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML247.M1.ID.NA[,1]) ~ NA),
         ZmCML247.M2 = case_when(ID %in% pull(ZmCML247.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML247.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML247.M2.ID.NA[,1]) ~ NA),
         ZmCML277.M1 = case_when(ID %in% pull(ZmCML277.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML277.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML277.M1.ID.NA[,1]) ~ NA),
         ZmCML277.M2 = case_when(ID %in% pull(ZmCML277.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML277.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML277.M2.ID.NA[,1]) ~ NA),
         ZmCML322.M1 = case_when(ID %in% pull(ZmCML322.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML322.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML322.M1.ID.NA[,1]) ~ NA),
         ZmCML322.M2 = case_when(ID %in% pull(ZmCML322.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML322.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML322.M2.ID.NA[,1]) ~ NA),
         ZmCML333.M1 = case_when(ID %in% pull(ZmCML333.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML333.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML333.M1.ID.NA[,1]) ~ NA),
         ZmCML333.M2 = case_when(ID %in% pull(ZmCML333.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML333.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML333.M2.ID.NA[,1]) ~ NA),
         ZmCML52.M1 = case_when(ID %in% pull(ZmCML52.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML52.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML52.M1.ID.NA[,1]) ~ NA),
         ZmCML52.M2 = case_when(ID %in% pull(ZmCML52.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML52.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML52.M2.ID.NA[,1]) ~ NA),
         ZmCML69.M1 = case_when(ID %in% pull(ZmCML69.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmCML69.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML69.M1.ID.NA[,1]) ~ NA),
         ZmCML69.M2 = case_when(ID %in% pull(ZmCML69.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmCML69.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmCML69.M2.ID.NA[,1]) ~ NA),
         ZmHP301.M1 = case_when(ID %in% pull(ZmHP301.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmHP301.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmHP301.M1.ID.NA[,1]) ~ NA),
         ZmHP301.M2 = case_when(ID %in% pull(ZmHP301.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmHP301.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmHP301.M2.ID.NA[,1]) ~ NA),
         ZmIL14H.M1 = case_when(ID %in% pull(ZmIL14H.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmIL14H.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmIL14H.M1.ID.NA[,1]) ~ NA),
         ZmIL14H.M2 = case_when(ID %in% pull(ZmIL14H.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmIL14H.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmIL14H.M2.ID.NA[,1]) ~ NA),
         ZmKi11.M1 = case_when(ID %in% pull(ZmKi11.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmKi11.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmKi11.M1.ID.NA[,1]) ~ NA),
         ZmKi11.M2 = case_when(ID %in% pull(ZmKi11.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmKi11.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmKi11.M2.ID.NA[,1]) ~ NA),
         ZmKi3.M1 = case_when(ID %in% pull(ZmKi3.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmKi3.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmKi3.M1.ID.NA[,1]) ~ NA),
         ZmKi3.M2 = case_when(ID %in% pull(ZmKi3.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmKi3.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmKi3.M2.ID.NA[,1]) ~ NA),
         ZmKy21.M1 = case_when(ID %in% pull(ZmKy21.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmKy21.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmKy21.M1.ID.NA[,1]) ~ NA),
         ZmKy21.M2 = case_when(ID %in% pull(ZmKy21.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmKy21.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmKy21.M2.ID.NA[,1]) ~ NA),
         ZmM162W.M1 = case_when(ID %in% pull(ZmM162W.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmM162W.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmM162W.M1.ID.NA[,1]) ~ NA),
         ZmM162W.M2 = case_when(ID %in% pull(ZmM162W.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmM162W.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmM162W.M2.ID.NA[,1]) ~ NA),
         ZmM37W.M1 = case_when(ID %in% pull(ZmM37W.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmM37W.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmM37W.M1.ID.NA[,1]) ~ NA),
         ZmM37W.M2 = case_when(ID %in% pull(ZmM37W.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmM37W.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmM37W.M2.ID.NA[,1]) ~ NA),
         ZmMo18W.M1 = case_when(ID %in% pull(ZmMo18W.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmMo18W.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmMo18W.M1.ID.NA[,1]) ~ NA),
         ZmMo18W.M2 = case_when(ID %in% pull(ZmMo18W.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmMo18W.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmMo18W.M2.ID.NA[,1]) ~ NA),
         ZmMS71.M1 = case_when(ID %in% pull(ZmMS71.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmMS71.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmMS71.M1.ID.NA[,1]) ~ NA),
         ZmMS71.M2 = case_when(ID %in% pull(ZmMS71.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmMS71.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmMS71.M2.ID.NA[,1]) ~ NA),
         ZmNC350.M1 = case_when(ID %in% pull(ZmNC350.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmNC350.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmNC350.M1.ID.NA[,1]) ~ NA),
         ZmNC350.M2 = case_when(ID %in% pull(ZmNC350.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmNC350.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmNC350.M2.ID.NA[,1]) ~ NA),
         ZmNC358.M1 = case_when(ID %in% pull(ZmNC358.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmNC358.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmNC358.M1.ID.NA[,1]) ~ NA),
         ZmNC358.M2 = case_when(ID %in% pull(ZmNC358.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmNC358.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmNC358.M2.ID.NA[,1]) ~ NA),
         ZmOh43.M1 = case_when(ID %in% pull(ZmOh43.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmOh43.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmOh43.M1.ID.NA[,1]) ~ NA),
         ZmOh43.M2 = case_when(ID %in% pull(ZmOh43.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmOh43.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmOh43.M2.ID.NA[,1]) ~ NA),
         ZmOh7b.M1 = case_when(ID %in% pull(ZmOh7b.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmOh7b.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmOh7b.M1.ID.NA[,1]) ~ NA),
         ZmOh7b.M2 = case_when(ID %in% pull(ZmOh7b.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmOh7b.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmOh7b.M2.ID.NA[,1]) ~ NA),
         ZmP39.M1 = case_when(ID %in% pull(ZmP39.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmP39.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmP39.M1.ID.NA[,1]) ~ NA),
         ZmP39.M2 = case_when(ID %in% pull(ZmP39.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmP39.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmP39.M2.ID.NA[,1]) ~ NA),
         ZmTx303.M1 = case_when(ID %in% pull(ZmTx303.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmTx303.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmTx303.M1.ID.NA[,1]) ~ NA),
         ZmTx303.M2 = case_when(ID %in% pull(ZmTx303.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmTx303.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmTx303.M2.ID.NA[,1]) ~ NA),
         ZmTzi8.M1 = case_when(ID %in% pull(ZmTzi8.M1.ID.del[,1]) ~ 1, ID %in% pull(ZmTzi8.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZmTzi8.M1.ID.NA[,1]) ~ NA),
         ZmTzi8.M2 = case_when(ID %in% pull(ZmTzi8.M2.ID.del[,1]) ~ 1, ID %in% pull(ZmTzi8.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZmTzi8.M2.ID.NA[,1]) ~ NA),
         ZnPI615697.M1 = case_when(ID %in% pull(ZnPI615697.M1.ID.del[,1]) ~ 1, ID %in% pull(ZnPI615697.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZnPI615697.M1.ID.NA[,1]) ~ NA),
         ZnPI615697.M2 = case_when(ID %in% pull(ZnPI615697.M2.ID.del[,1]) ~ 1, ID %in% pull(ZnPI615697.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZnPI615697.M2.ID.NA[,1]) ~ NA),
         ZvTIL01.M1 = case_when(ID %in% pull(ZvTIL01.M1.ID.del[,1]) ~ 1, ID %in% pull(ZvTIL01.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZvTIL01.M1.ID.NA[,1]) ~ NA),
         ZvTIL01.M2 = case_when(ID %in% pull(ZvTIL01.M2.ID.del[,1]) ~ 1, ID %in% pull(ZvTIL01.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZvTIL01.M2.ID.NA[,1]) ~ NA),
         ZvTIL11.M1 = case_when(ID %in% pull(ZvTIL11.M1.ID.del[,1]) ~ 1, ID %in% pull(ZvTIL11.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZvTIL11.M1.ID.NA[,1]) ~ NA),
         ZvTIL11.M2 = case_when(ID %in% pull(ZvTIL11.M2.ID.del[,1]) ~ 1, ID %in% pull(ZvTIL11.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZvTIL11.M2.ID.NA[,1]) ~ NA),
         ZxTIL18.M1 = case_when(ID %in% pull(ZxTIL18.M1.ID.del[,1]) ~ 1, ID %in% pull(ZxTIL18.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZxTIL18.M1.ID.NA[,1]) ~ NA),
         ZxTIL18.M2 = case_when(ID %in% pull(ZxTIL18.M2.ID.del[,1]) ~ 1, ID %in% pull(ZxTIL18.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZxTIL18.M2.ID.NA[,1]) ~ NA),
         ZxTIL25.M1 = case_when(ID %in% pull(ZxTIL25.M1.ID.del[,1]) ~ 1, ID %in% pull(ZxTIL25.M1.ID.ret[,1]) ~ 0, ID %in% pull(ZxTIL25.M1.ID.NA[,1]) ~ NA),
         ZxTIL25.M2 = case_when(ID %in% pull(ZxTIL25.M2.ID.del[,1]) ~ 1, ID %in% pull(ZxTIL25.M2.ID.ret[,1]) ~ 0, ID %in% pull(ZxTIL25.M2.ID.NA[,1]) ~ NA)
  )
#gives the total called fractionated across all Ms
summarize_if(full.fractionation.status, is.numeric,sum,na.rm=T) 

summarize_if(full.fractionation.status, is.numeric,sum,na.rm=T) %>% 
  pivot_longer(cols=colnames(.) ,names_to = "Subgenome", values_to = "Total_CDS_Fractionated")%>%
  mutate(genome = str_split(Subgenome, "[.]",simplify=T)[,1],
         M = str_split(Subgenome, "[.]",simplify=T)[,2],
         Percentage = (Total_CDS_Fractionated/69269)*100) %>%
  ggplot(aes(x=genome, y=Percentage))+
  geom_bar(aes(fill=M),stat="identity",position=position_dodge())+
  coord_flip()

#Check the overlap of calls vs. the genes called fractionated in NAM by Maggie
#NAM_standard<-tibble(M = c(rep(1,69),rep(2,57)),
#                     Gene_ID = c("Sobic.003G281800","Sobic.001G050400","Sobic.001G075700",
#                                 "Sobic.006G259201","Sobic.001G090700","Sobic.001G420800",
#                                 "Sobic.004G333600","Sobic.004G057100","Sobic.009G042800",
#                                 "Sobic.004G056400","Sobic.001G067900","Sobic.001G048832",
#                                 "Sobic.001G161000","Sobic.010G125000","Sobic.001G083800",
#                                 "Sobic.010G100800","Sobic.009G249500","Sobic.001G031300",
#                                 "Sobic.007G069000","Sobic.003G136200","Sobic.001G443300",
#                                 "Sobic.006G112100","Sobic.001G396200","Sobic.001G238700",
#                                 "Sobic.003G246000","Sobic.001G059900","Sobic.002G237700",
#                                 "Sobic.007G166800","Sobic.009G190500","Sobic.002G176200",
#                                 "Sobic.004G318200","Sobic.004G309600","Sobic.003G074300",
#                                 "Sobic.003G379000","Sobic.001G308900","Sobic.001G120500",
#                                 "Sobic.007G144600","Sobic.006G122900","Sobic.001G011700",
#                                 "Sobic.001G453800","Sobic.010G047700","Sobic.001G412400",
#                                 "Sobic.002G049300","Sobic.001G290700","Sobic.010G184000",
#                                 "Sobic.006G173900","Sobic.010G007100","Sobic.002G385300",
#                                 "Sobic.003G375200","Sobic.006G181600","Sobic.010G172400",
#                                 "Sobic.001G185800","Sobic.001G393300","Sobic.010G099700",
#                                 "Sobic.002G306400","Sobic.001G460300","Sobic.004G184200",
#                                 "Sobic.004G275200","Sobic.001G297900","Sobic.010G218300",
#                                 "Sobic.004G349000","Sobic.001G367700","Sobic.004G306200",
#                                 "Sobic.006G245400","Sobic.004G069200","Sobic.006G192100",
#                                 "Sobic.005G038700","Sobic.003G410400","Sobic.006G271000", #this is the break between M1 and M2 IDs
#                                 "Sobic.003G314100","Sobic.002G374700","Sobic.004G321800",
#                                 "Sobic.006G111000","Sobic.004G018100","Sobic.008G057000",
#                                 "Sobic.009G188000","Sobic.005G021000","Sobic.004G325100",
#                                "Sobic.002G257750","Sobic.002G017900","Sobic.010G023600",
#                                 "Sobic.006G201800","Sobic.009G210200","Sobic.009G024800",
#                                 "Sobic.003G253000","Sobic.002G139600","Sobic.006G200901",
#                                 "Sobic.009G221500","Sobic.006G197300","Sobic.010G015500",
#                                 "Sobic.001G295400","Sobic.010G224200","Sobic.006G085700",
#                                 "Sobic.001G410100","Sobic.003G389000","Sobic.002G214600",
#                                 "Sobic.009G216800","Sobic.006G203300","Sobic.010G044200",
#                                 "Sobic.004G305800","Sobic.003G392000","Sobic.006G121400",
#                                 "Sobic.010G088100","Sobic.004G286500","Sobic.001G308000",
 #                                "Sobic.009G006300","Sobic.010G218400","Sobic.003G253400",
#                                 "Sobic.006G095400","Sobic.004G321400","Sobic.003G074200",
#                                 "Sobic.010G236500","Sobic.008G018000","Sobic.001G290000",
#                                 "Sobic.002G371100","Sobic.009G149700","Sobic.003G139700",
#                                 "Sobic.003G214400","Sobic.001G365700","Sobic.006G082600",
#                                 "Sobic.004G211833","Sobic.006G258600","Sobic.006G112600",
#                                 "Sobic.009G061000","Sobic.006G269600","Sobic.003G221700"))
#NAM_standard.M1<-filter(NAM_standard, M == 1) 
#NAM_standard.M2<-filter(NAM_standard, M == 2)

#gene_fractionation %>% select(starts_with("Zm"), "Gene_ID") %>% select(ends_with(".M1"),"Gene_ID") %>%
#  filter(str_detect(Gene_ID, pattern = NAM_standard.M1$Gene_ID[2]))

#filter(ref_Sb313.cds, str_detect(Gene_ID, NAM_standard$Gene_ID[2]))

#agreement<-c()
#disagreement.1<-c()
#cantSay<-c()
#temp<-gene_fractionation %>% select(starts_with("Zm"), "Gene_ID") %>% select(ends_with(".M1"),"Gene_ID")
#for(i in 1:nrow(NAM_standard.M1)){
#  df<-filter(temp,str_detect(Gene_ID, pattern = NAM_standard.M1$Gene_ID[i]))
#  if(nrow(df)< 1){
#    cantSay<-c(cantSay,NAM_standard.M1$Gene_ID[i])
#  } else{
#    if(df %>% select(-Gene_ID) %>% rowSums() == 26){
#      agreement<-c(agreement,NAM_standard.M1$Gene_ID[i])
#    }else{
#      disagreement.1<-c(disagreement.1,NAM_standard.M1$Gene_ID[i])
#    }
#  }
#}
##M1 disagreement = 14, agreement = 33, cantSay (not in our ref set) = 22

#agreement<-c()
#disagreement<-c()
#cantSay<-c()
#temp<-gene_fractionation %>% select(starts_with("Zm"), "Gene_ID") %>% select(ends_with(".M2"),"Gene_ID")
#for(i in 1:nrow(NAM_standard.M2)){
#  df<-filter(temp,str_detect(Gene_ID, pattern = NAM_standard.M2$Gene_ID[i]))
#  if(nrow(df)< 1){
#    cantSay<-c(cantSay,NAM_standard.M2$Gene_ID[i])
#  } else{
#    if(df %>% select(-Gene_ID) %>% rowSums() == 26){
#      agreement<-c(agreement,NAM_standard.M2$Gene_ID[i])
#    }else{
#      disagreement<-c(disagreement,NAM_standard.M2$Gene_ID[i])
#    }
#  }
#}
##M2 disagreement = 14, agreement = 21, cantSay (not in our ref set) = 22

#NAM_standard<-read_csv("/work/LAS/mhufford-lab/snodgras/Fractionation/NAM_Standards.csv")
#any(is.na(NAM_standard$M2_Gene_ID))

#agreement<-c()
#disagreement<-c()
#cantSay<-c()
#temp<-gene_fractionation %>% select(starts_with("Zm"), "Gene_ID") %>% select(ends_with(".M1"),"Gene_ID")
#for(i in 1:nrow(NAM_standard)){
#  df<-filter(temp,str_detect(Gene_ID, pattern = NAM_standard$M1_Gene_ID[i]))
#  if(nrow(df)< 1){
#    cantSay<-c(cantSay,NAM_standard$M1_Gene_ID[i])
#  } else{
#    if(df %>% select(-Gene_ID) %>% rowSums() == 26){
#      agreement<-c(agreement,NAM_standard$M1_Gene_ID[i])
#    }else{
#      disagreement<-c(disagreement,NAM_standard$M1_Gene_ID[i])
#    }
#  }
#}
#2150 genes total
#agreement = 807, disagreement = 372, cantSay= 971

#agreement.2<-c()
#disagreement.2<-c()
#cantSay.2<-c()
#temp<-gene_fractionation %>% select(starts_with("Zm"), "Gene_ID") %>% select(ends_with(".M2"),"Gene_ID")
#for(i in 1:nrow(NAM_standard)){
#  if(!is.na(NAM_standard$M2_Gene_ID[i])){
#    df<-filter(temp,str_detect(Gene_ID, pattern = NAM_standard$M2_Gene_ID[i]))
#    if(nrow(df)< 1){
#      cantSay.2<-c(cantSay.2,NAM_standard$M2_Gene_ID[i])
#    } else{
#      if(df %>% select(-Gene_ID) %>% rowSums() == 26){
#        agreement.2<-c(agreement.2,NAM_standard$M2_Gene_ID[i])
#      }else{
#        disagreement.2<-c(disagreement.2,NAM_standard$M2_Gene_ID[i])
#      }
#    }
#  }
#}
#1701 genes total
#agreement = 668, disagreement = 298 , cantSay= 735 

#Categorical Fractionation Status

full.fractionation.status <-full.fractionation.status %>% mutate(TdFL.status = case_when(TdFL.M1 == 0 & TdFL.M2 == 0 ~ "Both_Retained",
                                                                                         TdFL.M1 == 1 & TdFL.M2 == 0 ~ "M2_Retained",
                                                                                         TdFL.M1 == 0 & TdFL.M2 == 1 ~ "M1_Retained",
                                                                                         TdFL.M1 == 1 & TdFL.M2 == 1 ~ "Both_Lost",
                                                                                         TdFL.M1 == 1 & is.na(TdFL.M2) ~ "M1_Lost:M2_NA",
                                                                                         TdFL.M1 == 0 & is.na(TdFL.M2) ~ "M1_Retained:M2_NA",
                                                                                         is.na(TdFL.M1) & TdFL.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                         is.na(TdFL.M1) & TdFL.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                         is.na(TdFL.M1) & is.na(TdFL.M2) ~ "Both_NA"),
                                                                 
                                                                 TdKS.status = case_when(TdKS.M1 == 0 & TdKS.M2 == 0 ~ "Both_Retained",
                                                                                         TdKS.M1 == 1 & TdKS.M2 == 0 ~ "M2_Retained",
                                                                                         TdKS.M1 == 0 & TdKS.M2 == 1 ~ "M1_Retained",
                                                                                         TdKS.M1 == 1 & TdKS.M2 == 1 ~ "Both_Lost",
                                                                                         TdKS.M1 == 1 & is.na(TdKS.M2) ~ "M1_Lost:M2_NA",
                                                                                         TdKS.M1 == 0 & is.na(TdKS.M2) ~ "M1_Retained:M2_NA",
                                                                                         is.na(TdKS.M1) & TdKS.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                         is.na(TdKS.M1) & TdKS.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                         is.na(TdKS.M1) & is.na(TdKS.M2) ~ "Both_NA"),
                                                                 
                                                                 ZdGigi_4to1.status = case_when(ZdGigi_4to1.M1 == 0 & ZdGigi_4to1.M2 == 0 ~ "Both_Retained",
                                                                                                ZdGigi_4to1.M1 == 1 & ZdGigi_4to1.M2 == 0 ~ "M2_Retained",
                                                                                                ZdGigi_4to1.M1 == 0 & ZdGigi_4to1.M2 == 1 ~ "M1_Retained",
                                                                                                ZdGigi_4to1.M1 == 1 & ZdGigi_4to1.M2 == 1 ~ "Both_Lost",
                                                                                                ZdGigi_4to1.M1 == 1 & is.na(ZdGigi_4to1.M2) ~ "M1_Lost:M2_NA",
                                                                                                ZdGigi_4to1.M1 == 0 & is.na(ZdGigi_4to1.M2) ~ "M1_Retained:M2_NA",
                                                                                         is.na(ZdGigi_4to1.M1) & ZdGigi_4to1.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                         is.na(ZdGigi_4to1.M1) & ZdGigi_4to1.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                         is.na(ZdGigi_4to1.M1) & is.na(ZdGigi_4to1.M2) ~ "Both_NA"),
                                                                 
                                                                 ZdMomo_4to1.status = case_when(ZdMomo_4to1.M1 == 0 & ZdMomo_4to1.M2 == 0 ~ "Both_Retained",
                                                                                                ZdMomo_4to1.M1 == 1 & ZdMomo_4to1.M2 == 0 ~ "M2_Retained",
                                                                                                ZdMomo_4to1.M1 == 0 & ZdMomo_4to1.M2 == 1 ~ "M1_Retained",
                                                                                                ZdMomo_4to1.M1 == 1 & ZdMomo_4to1.M2 == 1 ~ "Both_Lost",
                                                                                                ZdMomo_4to1.M1 == 1 & is.na(ZdMomo_4to1.M2) ~ "M1_Lost:M2_NA",
                                                                                                ZdMomo_4to1.M1 == 0 & is.na(ZdMomo_4to1.M2) ~ "M1_Retained:M2_NA",
                                                                                         is.na(ZdMomo_4to1.M1) & ZdMomo_4to1.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                         is.na(ZdMomo_4to1.M1) & ZdMomo_4to1.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                         is.na(ZdMomo_4to1.M1) & is.na(ZdMomo_4to1.M2) ~ "Both_NA"),
                                                                 
                                                                 ZnPI615697_4to1.status = case_when(ZnPI615697_4to1.M1 == 0 & ZnPI615697_4to1.M2 == 0 ~ "Both_Retained",
                                                                                                    ZnPI615697_4to1.M1 == 1 & ZnPI615697_4to1.M2 == 0 ~ "M2_Retained",
                                                                                                    ZnPI615697_4to1.M1 == 0 & ZnPI615697_4to1.M2 == 1 ~ "M1_Retained",
                                                                                                    ZnPI615697_4to1.M1 == 1 & ZnPI615697_4to1.M2 == 1 ~ "Both_Lost",
                                                                                                    ZnPI615697_4to1.M1 == 1 & is.na(ZnPI615697_4to1.M2) ~ "M1_Lost:M2_NA",
                                                                                                    ZnPI615697_4to1.M1 == 0 & is.na(ZnPI615697_4to1.M2) ~ "M1_Retained:M2_NA",
                                                                                         is.na(ZnPI615697_4to1.M1) & ZnPI615697_4to1.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                         is.na(ZnPI615697_4to1.M1) & ZnPI615697_4to1.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                         is.na(ZnPI615697_4to1.M1) & is.na(ZnPI615697_4to1.M2) ~ "Both_NA"),
                                                                 
                                                                 ZdGigi.status = case_when(ZdGigi.M1 == 0 & ZdGigi.M2 == 0 ~ "Both_Retained",
                                                                                           ZdGigi.M1 == 1 & ZdGigi.M2 == 0 ~ "M2_Retained",
                                                                                           ZdGigi.M1 == 0 & ZdGigi.M2 == 1 ~ "M1_Retained",
                                                                                           ZdGigi.M1 == 1 & ZdGigi.M2 == 1 ~ "Both_Lost",
                                                                                           ZdGigi.M1 == 1 & is.na(ZdGigi.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZdGigi.M1 == 0 & is.na(ZdGigi.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZdGigi.M1) & ZdGigi.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZdGigi.M1) & ZdGigi.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZdGigi.M1) & is.na(ZdGigi.M2) ~ "Both_NA"),
                                                                 
                                                                 ZdMomo.status = case_when(ZdMomo.M1 == 0 & ZdMomo.M2 == 0 ~ "Both_Retained",
                                                                                           ZdMomo.M1 == 1 & ZdMomo.M2 == 0 ~ "M2_Retained",
                                                                                           ZdMomo.M1 == 0 & ZdMomo.M2 == 1 ~ "M1_Retained",
                                                                                           ZdMomo.M1 == 1 & ZdMomo.M2 == 1 ~ "Both_Lost",
                                                                                           ZdMomo.M1 == 1 & is.na(ZdMomo.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZdMomo.M1 == 0 & is.na(ZdMomo.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZdMomo.M1) & ZdMomo.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZdMomo.M1) & ZdMomo.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZdMomo.M1) & is.na(ZdMomo.M2) ~ "Both_NA"),
                                                                 
                                                                 ZhRIMHU001.status = case_when(ZhRIMHU001.M1 == 0 & ZhRIMHU001.M2 == 0 ~ "Both_Retained",
                                                                                               ZhRIMHU001.M1 == 1 & ZhRIMHU001.M2 == 0 ~ "M2_Retained",
                                                                                               ZhRIMHU001.M1 == 0 & ZhRIMHU001.M2 == 1 ~ "M1_Retained",
                                                                                               ZhRIMHU001.M1 == 1 & ZhRIMHU001.M2 == 1 ~ "Both_Lost",
                                                                                               ZhRIMHU001.M1 == 1 & is.na(ZhRIMHU001.M2) ~ "M1_Lost:M2_NA",
                                                                                               ZhRIMHU001.M1 == 0 & is.na(ZhRIMHU001.M2) ~ "M1_Retained:M2_NA",
                                                                                               is.na(ZhRIMHU001.M1) & ZhRIMHU001.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                               is.na(ZhRIMHU001.M1) & ZhRIMHU001.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                               is.na(ZhRIMHU001.M1) & is.na(ZhRIMHU001.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmB73.status = case_when(ZmB73.M1 == 0 & ZmB73.M2 == 0 ~ "Both_Retained",
                                                                                          ZmB73.M1 == 1 & ZmB73.M2 == 0 ~ "M2_Retained",
                                                                                          ZmB73.M1 == 0 & ZmB73.M2 == 1 ~ "M1_Retained",
                                                                                          ZmB73.M1 == 1 & ZmB73.M2 == 1 ~ "Both_Lost",
                                                                                          ZmB73.M1 == 1 & is.na(ZmB73.M2) ~ "M1_Lost:M2_NA",
                                                                                          ZmB73.M1 == 0 & is.na(ZmB73.M2) ~ "M1_Retained:M2_NA",
                                                                                          is.na(ZmB73.M1) & ZmB73.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                          is.na(ZmB73.M1) & ZmB73.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                          is.na(ZmB73.M1) & is.na(ZmB73.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmB97.status = case_when(ZmB97.M1 == 0 & ZmB97.M2 == 0 ~ "Both_Retained",
                                                                                          ZmB97.M1 == 1 & ZmB97.M2 == 0 ~ "M2_Retained",
                                                                                          ZmB97.M1 == 0 & ZmB97.M2 == 1 ~ "M1_Retained",
                                                                                          ZmB97.M1 == 1 & ZmB97.M2 == 1 ~ "Both_Lost",
                                                                                          ZmB97.M1 == 1 & is.na(ZmB97.M2) ~ "M1_Lost:M2_NA",
                                                                                          ZmB97.M1 == 0 & is.na(ZmB97.M2) ~ "M1_Retained:M2_NA",
                                                                                          is.na(ZmB97.M1) & ZmB97.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                          is.na(ZmB97.M1) & ZmB97.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                          is.na(ZmB97.M1) & is.na(ZmB97.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML103.status = case_when(ZmCML103.M1 == 0 & ZmCML103.M2 == 0 ~ "Both_Retained",
                                                                                             ZmCML103.M1 == 1 & ZmCML103.M2 == 0 ~ "M2_Retained",
                                                                                             ZmCML103.M1 == 0 & ZmCML103.M2 == 1 ~ "M1_Retained",
                                                                                             ZmCML103.M1 == 1 & ZmCML103.M2 == 1 ~ "Both_Lost",
                                                                                             ZmCML103.M1 == 1 & is.na(ZmCML103.M2) ~ "M1_Lost:M2_NA",
                                                                                             ZmCML103.M1 == 0 & is.na(ZmCML103.M2) ~ "M1_Retained:M2_NA",
                                                                                             is.na(ZmCML103.M1) & ZmCML103.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                             is.na(ZmCML103.M1) & ZmCML103.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                             is.na(ZmCML103.M1) & is.na(ZmCML103.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML228.status = case_when(ZmCML228.M1 == 0 & ZmCML228.M2 == 0 ~ "Both_Retained",
                                                                                             ZmCML228.M1 == 1 & ZmCML228.M2 == 0 ~ "M2_Retained",
                                                                                             ZmCML228.M1 == 0 & ZmCML228.M2 == 1 ~ "M1_Retained",
                                                                                             ZmCML228.M1 == 1 & ZmCML228.M2 == 1 ~ "Both_Lost",
                                                                                             ZmCML228.M1 == 1 & is.na(ZmCML228.M2) ~ "M1_Lost:M2_NA",
                                                                                             ZmCML228.M1 == 0 & is.na(ZmCML228.M2) ~ "M1_Retained:M2_NA",
                                                                                             is.na(ZmCML228.M1) & ZmCML228.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                             is.na(ZmCML228.M1) & ZmCML228.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                             is.na(ZmCML228.M1) & is.na(ZmCML228.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML247.status = case_when(ZmCML247.M1 == 0 & ZmCML247.M2 == 0 ~ "Both_Retained",
                                                                                             ZmCML247.M1 == 1 & ZmCML247.M2 == 0 ~ "M2_Retained",
                                                                                             ZmCML247.M1 == 0 & ZmCML247.M2 == 1 ~ "M1_Retained",
                                                                                             ZmCML247.M1 == 1 & ZmCML247.M2 == 1 ~ "Both_Lost",
                                                                                             ZmCML247.M1 == 1 & is.na(ZmCML247.M2) ~ "M1_Lost:M2_NA",
                                                                                             ZmCML247.M1 == 0 & is.na(ZmCML247.M2) ~ "M1_Retained:M2_NA",
                                                                                             is.na(ZmCML247.M1) & ZmCML247.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                             is.na(ZmCML247.M1) & ZmCML247.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                             is.na(ZmCML247.M1) & is.na(ZmCML247.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML277.status = case_when(ZmCML277.M1 == 0 & ZmCML277.M2 == 0 ~ "Both_Retained",
                                                                                             ZmCML277.M1 == 1 & ZmCML277.M2 == 0 ~ "M2_Retained",
                                                                                             ZmCML277.M1 == 0 & ZmCML277.M2 == 1 ~ "M1_Retained",
                                                                                             ZmCML277.M1 == 1 & ZmCML277.M2 == 1 ~ "Both_Lost",
                                                                                             ZmCML277.M1 == 1 & is.na(ZmCML277.M2) ~ "M1_Lost:M2_NA",
                                                                                             ZmCML277.M1 == 0 & is.na(ZmCML277.M2) ~ "M1_Retained:M2_NA",
                                                                                             is.na(ZmCML277.M1) & ZmCML277.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                             is.na(ZmCML277.M1) & ZmCML277.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                             is.na(ZmCML277.M1) & is.na(ZmCML277.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML322.status = case_when(ZmCML322.M1 == 0 & ZmCML322.M2 == 0 ~ "Both_Retained",
                                                                                             ZmCML322.M1 == 1 & ZmCML322.M2 == 0 ~ "M2_Retained",
                                                                                             ZmCML322.M1 == 0 & ZmCML322.M2 == 1 ~ "M1_Retained",
                                                                                             ZmCML322.M1 == 1 & ZmCML322.M2 == 1 ~ "Both_Lost",
                                                                                             ZmCML322.M1 == 1 & is.na(ZmCML322.M2) ~ "M1_Lost:M2_NA",
                                                                                             ZmCML322.M1 == 0 & is.na(ZmCML322.M2) ~ "M1_Retained:M2_NA",
                                                                                             is.na(ZmCML322.M1) & ZmCML322.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                             is.na(ZmCML322.M1) & ZmCML322.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                             is.na(ZmCML322.M1) & is.na(ZmCML322.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML333.status = case_when(ZmCML333.M1 == 0 & ZmCML333.M2 == 0 ~ "Both_Retained",
                                                                                             ZmCML333.M1 == 1 & ZmCML333.M2 == 0 ~ "M2_Retained",
                                                                                             ZmCML333.M1 == 0 & ZmCML333.M2 == 1 ~ "M1_Retained",
                                                                                             ZmCML333.M1 == 1 & ZmCML333.M2 == 1 ~ "Both_Lost",
                                                                                             ZmCML333.M1 == 1 & is.na(ZmCML333.M2) ~ "M1_Lost:M2_NA",
                                                                                             ZmCML333.M1 == 0 & is.na(ZmCML333.M2) ~ "M1_Retained:M2_NA",
                                                                                             is.na(ZmCML333.M1) & ZmCML333.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                             is.na(ZmCML333.M1) & ZmCML333.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                             is.na(ZmCML333.M1) & is.na(ZmCML333.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML52.status = case_when(ZmCML52.M1 == 0 & ZmCML52.M2 == 0 ~ "Both_Retained",
                                                                                            ZmCML52.M1 == 1 & ZmCML52.M2 == 0 ~ "M2_Retained",
                                                                                            ZmCML52.M1 == 0 & ZmCML52.M2 == 1 ~ "M1_Retained",
                                                                                            ZmCML52.M1 == 1 & ZmCML52.M2 == 1 ~ "Both_Lost",
                                                                                            ZmCML52.M1 == 1 & is.na(ZmCML52.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmCML52.M1 == 0 & is.na(ZmCML52.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmCML52.M1) & ZmCML52.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmCML52.M1) & ZmCML52.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmCML52.M1) & is.na(ZmCML52.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmCML69.status = case_when(ZmCML69.M1 == 0 & ZmCML69.M2 == 0 ~ "Both_Retained",
                                                                                            ZmCML69.M1 == 1 & ZmCML69.M2 == 0 ~ "M2_Retained",
                                                                                            ZmCML69.M1 == 0 & ZmCML69.M2 == 1 ~ "M1_Retained",
                                                                                            ZmCML69.M1 == 1 & ZmCML69.M2 == 1 ~ "Both_Lost",
                                                                                            ZmCML69.M1 == 1 & is.na(ZmCML69.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmCML69.M1 == 0 & is.na(ZmCML69.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmCML69.M1) & ZmCML69.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmCML69.M1) & ZmCML69.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmCML69.M1) & is.na(ZmCML69.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmHP301.status = case_when(ZmHP301.M1 == 0 & ZmHP301.M2 == 0 ~ "Both_Retained",
                                                                                            ZmHP301.M1 == 1 & ZmHP301.M2 == 0 ~ "M2_Retained",
                                                                                            ZmHP301.M1 == 0 & ZmHP301.M2 == 1 ~ "M1_Retained",
                                                                                            ZmHP301.M1 == 1 & ZmHP301.M2 == 1 ~ "Both_Lost",
                                                                                            ZmHP301.M1 == 1 & is.na(ZmHP301.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmHP301.M1 == 0 & is.na(ZmHP301.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmHP301.M1) & ZmHP301.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmHP301.M1) & ZmHP301.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmHP301.M1) & is.na(ZmHP301.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmIL14H.status = case_when(ZmIL14H.M1 == 0 & ZmIL14H.M2 == 0 ~ "Both_Retained",
                                                                                            ZmIL14H.M1 == 1 & ZmIL14H.M2 == 0 ~ "M2_Retained",
                                                                                            ZmIL14H.M1 == 0 & ZmIL14H.M2 == 1 ~ "M1_Retained",
                                                                                            ZmIL14H.M1 == 1 & ZmIL14H.M2 == 1 ~ "Both_Lost",
                                                                                            ZmIL14H.M1 == 1 & is.na(ZmIL14H.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmIL14H.M1 == 0 & is.na(ZmIL14H.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmIL14H.M1) & ZmIL14H.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmIL14H.M1) & ZmIL14H.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmIL14H.M1) & is.na(ZmIL14H.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmKi11.status = case_when(ZmKi11.M1 == 0 & ZmKi11.M2 == 0 ~ "Both_Retained",
                                                                                           ZmKi11.M1 == 1 & ZmKi11.M2 == 0 ~ "M2_Retained",
                                                                                           ZmKi11.M1 == 0 & ZmKi11.M2 == 1 ~ "M1_Retained",
                                                                                           ZmKi11.M1 == 1 & ZmKi11.M2 == 1 ~ "Both_Lost",
                                                                                           ZmKi11.M1 == 1 & is.na(ZmKi11.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmKi11.M1 == 0 & is.na(ZmKi11.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmKi11.M1) & ZmKi11.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmKi11.M1) & ZmKi11.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmKi11.M1) & is.na(ZmKi11.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmKi3.status = case_when(ZmKi3.M1 == 0 & ZmKi3.M2 == 0 ~ "Both_Retained",
                                                                                          ZmKi3.M1 == 1 & ZmKi3.M2 == 0 ~ "M2_Retained",
                                                                                          ZmKi3.M1 == 0 & ZmKi3.M2 == 1 ~ "M1_Retained",
                                                                                          ZmKi3.M1 == 1 & ZmKi3.M2 == 1 ~ "Both_Lost",
                                                                                          ZmKi3.M1 == 1 & is.na(ZmKi3.M2) ~ "M1_Lost:M2_NA",
                                                                                          ZmKi3.M1 == 0 & is.na(ZmKi3.M2) ~ "M1_Retained:M2_NA",
                                                                                          is.na(ZmKi3.M1) & ZmKi3.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                          is.na(ZmKi3.M1) & ZmKi3.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                          is.na(ZmKi3.M1) & is.na(ZmKi3.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmKy21.status = case_when(ZmKy21.M1 == 0 & ZmKy21.M2 == 0 ~ "Both_Retained",
                                                                                           ZmKy21.M1 == 1 & ZmKy21.M2 == 0 ~ "M2_Retained",
                                                                                           ZmKy21.M1 == 0 & ZmKy21.M2 == 1 ~ "M1_Retained",
                                                                                           ZmKy21.M1 == 1 & ZmKy21.M2 == 1 ~ "Both_Lost",
                                                                                           ZmKy21.M1 == 1 & is.na(ZmKy21.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmKy21.M1 == 0 & is.na(ZmKy21.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmKy21.M1) & ZmKy21.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmKy21.M1) & ZmKy21.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmKy21.M1) & is.na(ZmKy21.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmM162W.status = case_when(ZmM162W.M1 == 0 & ZmM162W.M2 == 0 ~ "Both_Retained",
                                                                                            ZmM162W.M1 == 1 & ZmM162W.M2 == 0 ~ "M2_Retained",
                                                                                            ZmM162W.M1 == 0 & ZmM162W.M2 == 1 ~ "M1_Retained",
                                                                                            ZmM162W.M1 == 1 & ZmM162W.M2 == 1 ~ "Both_Lost",
                                                                                            ZmM162W.M1 == 1 & is.na(ZmM162W.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmM162W.M1 == 0 & is.na(ZmM162W.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmM162W.M1) & ZmM162W.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmM162W.M1) & ZmM162W.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmM162W.M1) & is.na(ZmM162W.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmM37W.status = case_when(ZmM37W.M1 == 0 & ZmM37W.M2 == 0 ~ "Both_Retained",
                                                                                           ZmM37W.M1 == 1 & ZmM37W.M2 == 0 ~ "M2_Retained",
                                                                                           ZmM37W.M1 == 0 & ZmM37W.M2 == 1 ~ "M1_Retained",
                                                                                           ZmM37W.M1 == 1 & ZmM37W.M2 == 1 ~ "Both_Lost",
                                                                                           ZmM37W.M1 == 1 & is.na(ZmM37W.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmM37W.M1 == 0 & is.na(ZmM37W.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmM37W.M1) & ZmM37W.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmM37W.M1) & ZmM37W.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmM37W.M1) & is.na(ZmM37W.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmMo18W.status = case_when(ZmMo18W.M1 == 0 & ZmMo18W.M2 == 0 ~ "Both_Retained",
                                                                                            ZmMo18W.M1 == 1 & ZmMo18W.M2 == 0 ~ "M2_Retained",
                                                                                            ZmMo18W.M1 == 0 & ZmMo18W.M2 == 1 ~ "M1_Retained",
                                                                                            ZmMo18W.M1 == 1 & ZmMo18W.M2 == 1 ~ "Both_Lost",
                                                                                            ZmMo18W.M1 == 1 & is.na(ZmMo18W.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmMo18W.M1 == 0 & is.na(ZmMo18W.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmMo18W.M1) & ZmMo18W.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmMo18W.M1) & ZmMo18W.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmMo18W.M1) & is.na(ZmMo18W.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmMS71.status = case_when(ZmMS71.M1 == 0 & ZmMS71.M2 == 0 ~ "Both_Retained",
                                                                                           ZmMS71.M1 == 1 & ZmMS71.M2 == 0 ~ "M2_Retained",
                                                                                           ZmMS71.M1 == 0 & ZmMS71.M2 == 1 ~ "M1_Retained",
                                                                                           ZmMS71.M1 == 1 & ZmMS71.M2 == 1 ~ "Both_Lost",
                                                                                           ZmMS71.M1 == 1 & is.na(ZmMS71.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmMS71.M1 == 0 & is.na(ZmMS71.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmMS71.M1) & ZmMS71.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmMS71.M1) & ZmMS71.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmMS71.M1) & is.na(ZmMS71.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmNC350.status = case_when(ZmNC350.M1 == 0 & ZmNC350.M2 == 0 ~ "Both_Retained",
                                                                                            ZmNC350.M1 == 1 & ZmNC350.M2 == 0 ~ "M2_Retained",
                                                                                            ZmNC350.M1 == 0 & ZmNC350.M2 == 1 ~ "M1_Retained",
                                                                                            ZmNC350.M1 == 1 & ZmNC350.M2 == 1 ~ "Both_Lost",
                                                                                            ZmNC350.M1 == 1 & is.na(ZmNC350.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmNC350.M1 == 0 & is.na(ZmNC350.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmNC350.M1) & ZmNC350.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmNC350.M1) & ZmNC350.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmNC350.M1) & is.na(ZmNC350.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmNC358.status = case_when(ZmNC358.M1 == 0 & ZmNC358.M2 == 0 ~ "Both_Retained",
                                                                                            ZmNC358.M1 == 1 & ZmNC358.M2 == 0 ~ "M2_Retained",
                                                                                            ZmNC358.M1 == 0 & ZmNC358.M2 == 1 ~ "M1_Retained",
                                                                                            ZmNC358.M1 == 1 & ZmNC358.M2 == 1 ~ "Both_Lost",
                                                                                            ZmNC358.M1 == 1 & is.na(ZmNC358.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmNC358.M1 == 0 & is.na(ZmNC358.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmNC358.M1) & ZmNC358.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmNC358.M1) & ZmNC358.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmNC358.M1) & is.na(ZmNC358.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmOh43.status = case_when(ZmOh43.M1 == 0 & ZmOh43.M2 == 0 ~ "Both_Retained",
                                                                                           ZmOh43.M1 == 1 & ZmOh43.M2 == 0 ~ "M2_Retained",
                                                                                           ZmOh43.M1 == 0 & ZmOh43.M2 == 1 ~ "M1_Retained",
                                                                                           ZmOh43.M1 == 1 & ZmOh43.M2 == 1 ~ "Both_Lost",
                                                                                           ZmOh43.M1 == 1 & is.na(ZmOh43.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmOh43.M1 == 0 & is.na(ZmOh43.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmOh43.M1) & ZmOh43.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmOh43.M1) & ZmOh43.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmOh43.M1) & is.na(ZmOh43.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmOh7b.status = case_when(ZmOh7b.M1 == 0 & ZmOh7b.M2 == 0 ~ "Both_Retained",
                                                                                           ZmOh7b.M1 == 1 & ZmOh7b.M2 == 0 ~ "M2_Retained",
                                                                                           ZmOh7b.M1 == 0 & ZmOh7b.M2 == 1 ~ "M1_Retained",
                                                                                           ZmOh7b.M1 == 1 & ZmOh7b.M2 == 1 ~ "Both_Lost",
                                                                                           ZmOh7b.M1 == 1 & is.na(ZmOh7b.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmOh7b.M1 == 0 & is.na(ZmOh7b.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmOh7b.M1) & ZmOh7b.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmOh7b.M1) & ZmOh7b.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmOh7b.M1) & is.na(ZmOh7b.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmP39.status = case_when(ZmP39.M1 == 0 & ZmP39.M2 == 0 ~ "Both_Retained",
                                                                                          ZmP39.M1 == 1 & ZmP39.M2 == 0 ~ "M2_Retained",
                                                                                          ZmP39.M1 == 0 & ZmP39.M2 == 1 ~ "M1_Retained",
                                                                                          ZmP39.M1 == 1 & ZmP39.M2 == 1 ~ "Both_Lost",
                                                                                          ZmP39.M1 == 1 & is.na(ZmP39.M2) ~ "M1_Lost:M2_NA",
                                                                                          ZmP39.M1 == 0 & is.na(ZmP39.M2) ~ "M1_Retained:M2_NA",
                                                                                          is.na(ZmP39.M1) & ZmP39.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                          is.na(ZmP39.M1) & ZmP39.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                          is.na(ZmP39.M1) & is.na(ZmP39.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmTx303.status = case_when(ZmTx303.M1 == 0 & ZmTx303.M2 == 0 ~ "Both_Retained",
                                                                                            ZmTx303.M1 == 1 & ZmTx303.M2 == 0 ~ "M2_Retained",
                                                                                            ZmTx303.M1 == 0 & ZmTx303.M2 == 1 ~ "M1_Retained",
                                                                                            ZmTx303.M1 == 1 & ZmTx303.M2 == 1 ~ "Both_Lost",
                                                                                            ZmTx303.M1 == 1 & is.na(ZmTx303.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZmTx303.M1 == 0 & is.na(ZmTx303.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZmTx303.M1) & ZmTx303.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZmTx303.M1) & ZmTx303.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZmTx303.M1) & is.na(ZmTx303.M2) ~ "Both_NA"),
                                                                 
                                                                 ZmTzi8.status = case_when(ZmTzi8.M1 == 0 & ZmTzi8.M2 == 0 ~ "Both_Retained",
                                                                                           ZmTzi8.M1 == 1 & ZmTzi8.M2 == 0 ~ "M2_Retained",
                                                                                           ZmTzi8.M1 == 0 & ZmTzi8.M2 == 1 ~ "M1_Retained",
                                                                                           ZmTzi8.M1 == 1 & ZmTzi8.M2 == 1 ~ "Both_Lost",
                                                                                           ZmTzi8.M1 == 1 & is.na(ZmTzi8.M2) ~ "M1_Lost:M2_NA",
                                                                                           ZmTzi8.M1 == 0 & is.na(ZmTzi8.M2) ~ "M1_Retained:M2_NA",
                                                                                           is.na(ZmTzi8.M1) & ZmTzi8.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                           is.na(ZmTzi8.M1) & ZmTzi8.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                           is.na(ZmTzi8.M1) & is.na(ZmTzi8.M2) ~ "Both_NA"),
                                                                 
                                                                 ZnPI615697.status = case_when(ZnPI615697.M1 == 0 & ZnPI615697.M2 == 0 ~ "Both_Retained",
                                                                                               ZnPI615697.M1 == 1 & ZnPI615697.M2 == 0 ~ "M2_Retained",
                                                                                               ZnPI615697.M1 == 0 & ZnPI615697.M2 == 1 ~ "M1_Retained",
                                                                                               ZnPI615697.M1 == 1 & ZnPI615697.M2 == 1 ~ "Both_Lost",
                                                                                               ZnPI615697.M1 == 1 & is.na(ZnPI615697.M2) ~ "M1_Lost:M2_NA",
                                                                                               ZnPI615697.M1 == 0 & is.na(ZnPI615697.M2) ~ "M1_Retained:M2_NA",
                                                                                               is.na(ZnPI615697.M1) & ZnPI615697.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                               is.na(ZnPI615697.M1) & ZnPI615697.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                               is.na(ZnPI615697.M1) & is.na(ZnPI615697.M2) ~ "Both_NA"),
                                                                 
                                                                 ZvTIL01.status = case_when(ZvTIL01.M1 == 0 & ZvTIL01.M2 == 0 ~ "Both_Retained",
                                                                                            ZvTIL01.M1 == 1 & ZvTIL01.M2 == 0 ~ "M2_Retained",
                                                                                            ZvTIL01.M1 == 0 & ZvTIL01.M2 == 1 ~ "M1_Retained",
                                                                                            ZvTIL01.M1 == 1 & ZvTIL01.M2 == 1 ~ "Both_Lost",
                                                                                            ZvTIL01.M1 == 1 & is.na(ZvTIL01.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZvTIL01.M1 == 0 & is.na(ZvTIL01.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZvTIL01.M1) & ZvTIL01.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZvTIL01.M1) & ZvTIL01.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZvTIL01.M1) & is.na(ZvTIL01.M2) ~ "Both_NA"),
                                                                 
                                                                 ZvTIL11.status = case_when(ZvTIL11.M1 == 0 & ZvTIL11.M2 == 0 ~ "Both_Retained",
                                                                                            ZvTIL11.M1 == 1 & ZvTIL11.M2 == 0 ~ "M2_Retained",
                                                                                            ZvTIL11.M1 == 0 & ZvTIL11.M2 == 1 ~ "M1_Retained",
                                                                                            ZvTIL11.M1 == 1 & ZvTIL11.M2 == 1 ~ "Both_Lost",
                                                                                            ZvTIL11.M1 == 1 & is.na(ZvTIL11.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZvTIL11.M1 == 0 & is.na(ZvTIL11.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZvTIL11.M1) & ZvTIL11.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZvTIL11.M1) & ZvTIL11.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZvTIL11.M1) & is.na(ZvTIL11.M2) ~ "Both_NA"),
                                                                 
                                                                 ZxTIL18.status = case_when(ZxTIL18.M1 == 0 & ZxTIL18.M2 == 0 ~ "Both_Retained",
                                                                                            ZxTIL18.M1 == 1 & ZxTIL18.M2 == 0 ~ "M2_Retained",
                                                                                            ZxTIL18.M1 == 0 & ZxTIL18.M2 == 1 ~ "M1_Retained",
                                                                                            ZxTIL18.M1 == 1 & ZxTIL18.M2 == 1 ~ "Both_Lost",
                                                                                            ZxTIL18.M1 == 1 & is.na(ZxTIL18.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZxTIL18.M1 == 0 & is.na(ZxTIL18.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZxTIL18.M1) & ZxTIL18.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZxTIL18.M1) & ZxTIL18.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZxTIL18.M1) & is.na(ZxTIL18.M2) ~ "Both_NA"),
                                                                 
                                                                 ZxTIL25.status = case_when(ZxTIL25.M1 == 0 & ZxTIL25.M2 == 0 ~ "Both_Retained",
                                                                                            ZxTIL25.M1 == 1 & ZxTIL25.M2 == 0 ~ "M2_Retained",
                                                                                            ZxTIL25.M1 == 0 & ZxTIL25.M2 == 1 ~ "M1_Retained",
                                                                                            ZxTIL25.M1 == 1 & ZxTIL25.M2 == 1 ~ "Both_Lost",
                                                                                            ZxTIL25.M1 == 1 & is.na(ZxTIL25.M2) ~ "M1_Lost:M2_NA",
                                                                                            ZxTIL25.M1 == 0 & is.na(ZxTIL25.M2) ~ "M1_Retained:M2_NA",
                                                                                            is.na(ZxTIL25.M1) & ZxTIL25.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                            is.na(ZxTIL25.M1) & ZxTIL25.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                            is.na(ZxTIL25.M1) & is.na(ZxTIL25.M2) ~ "Both_NA"),
                                                                 
)
long.full.fractionation.status <-full.fractionation.status %>% dplyr::select(c(ID,contains("status"))) %>% 
  pivot_longer(cols = ends_with("status"), names_to = "Genome", values_to = "Status" )
long.full.fractionation.status$Genome<-long.full.fractionation.status$Genome %>% str_remove(".status")
#Tripsacum = "#F45B69", #Zn ="#114B5F", #Zd = "#2B5A78" #Zh = "#247590" #Zmex = "#03A0B5" #Zparv = "#4DA89D" 
#Zmays TROPICAL = "#7CCBB2" MIXED+HP301 = "#92DDB0" TEMP+Sweet = "#DAFDD8"

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
                 "TdFL" = "#F45B69",
                 "ZdGigi_4to1" = "#433475",
                 "ZdMomo_4to1" = "#433475",
                 "ZnPI615697_4to1" = "#412151",
                 "TdKS" = "#F45B69")
long.full.fractionation.status$Genome<-long.full.fractionation.status$Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                                                                                 "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                                                                                 "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                                                                                 "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))
long.full.fractionation.status$Status<-long.full.fractionation.status$Status %>% factor(levels = c("Both_Retained","M1_Retained","M2_Retained", "Both_Lost","M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA"))
ggplot(data = long.full.fractionation.status, aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Number of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.AllGenomes.png", device="png",dpi=500,width = 8, height = 7)

ggplot(data = filter(long.full.fractionation.status, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697")), aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Number of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.ZnZd1to4only.png", device="png",dpi=500,width = 8, height = 7)

ggplot(data = filter(long.full.fractionation.status, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1")), aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Number of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.noZnZd.png", device="png",dpi=500,width = 8, height = 7)

group_by(long.full.fractionation.status, Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(y=Status, x=Percent))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Percentage of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.AllGenomes.png", device="png",dpi=500,width = 8, height = 7)

filter(long.full.fractionation.status, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697")) %>% group_by(Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(y=Status, x=Percent))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Percentage of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.ZnZd4to1only.png", device="png",dpi=500,width = 8, height = 7)

filter(long.full.fractionation.status, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1")) %>% group_by(Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(y=Status, x=Percent))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Percentage of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.NoZnZd.png", device="png",dpi=500,width = 8, height = 7)

filter(long.full.fractionation.status, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1") & Status %in% c("Both_Retained","M1_Retained","M2_Retained","Both_Lost")) %>% group_by(Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(y=Status, x=Percent))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Percentage of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Deleted","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.NoZnZd.AlignedOnly.png", device="png",dpi=500,width = 8, height = 7)

filter(long.full.fractionation.status, !Genome %in% c("ZdGigi","ZdMomo","ZnPI615697") & Status %in% c("Both_Retained","M1_Retained","M2_Retained","Both_Lost")) %>% group_by(Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(y=Status, x=Percent))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Percentage of Exons")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Deleted","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("")+
  ggtitle("Fractionation Status by Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.ZnZd4to1only.AlignedOnly.png", device="png",dpi=500,width = 8, height = 7)


filter(long.full.fractionation.status, Status == "Both_Retained") %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#13884  7174 19267
filter(long.full.fractionation.status, Status %in% c("M1_Retained","M2_Retained")) %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#31763. 16643 33636
filter(long.full.fractionation.status, Status == "Both_Lost") %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#14069.  7542 14973
filter(long.full.fractionation.status, Status %in% c("M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA")) %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#9552.  5952 37910

#Summary stats without Zn and only Zd 4to1
filter(long.full.fractionation.status, Status == "Both_Retained" & !Genome %in% c("ZnPI615697", "ZnPI615697_4to1","ZdMomo","ZdGigi")) %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#14342. 13290 19267
filter(long.full.fractionation.status, Status %in% c("M1_Retained","M2_Retained") & !Genome %in% c("ZnPI615697_4to1","ZnPI615697", "ZdMomo","ZdGigi")) %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#32804. 26988 33636
filter(long.full.fractionation.status, Status == "Both_Lost" & !Genome %in% c("ZnPI615697_4to1","ZnPI615697", "ZdMomo","ZdGigi")) %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#14541. 12497 14973
filter(long.full.fractionation.status, Status %in% c("M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA") & !Genome %in% c("ZnPI615697_4to1","ZnPI615697", "ZdMomo","ZdGigi")) %>% group_by(Genome) %>% count() %>% as_tibble() %>%
  reframe(avg.n = mean(n), min=min(n),max=max(n))
#7582.  5952  13049

##Summarize Fractionation Status by gene
full.fractionation.status<-full.fractionation.status %>% 
  mutate(CDS_ID = str_split(ID, ";", simplify=T)[,1] %>% str_remove_all(pattern = "ID="),
         Gene_ID = str_split(ID, ";",simplify = T)[,2] %>% str_remove_all(pattern = "Parent="))

gene_id<-select(full.fractionation.status, Gene_ID) %>% pull() %>% unique()

#Start making gene fractionation status by summing the statuses for each exon for a given gene
#sum will be # of exons that are called fractionated for a given gene
gene_fractionation<-tibble(Gene_ID=NA, TdFL.M1.sum=NA,TdFL.M2.sum=NA,ZdGigi.M1.sum=NA,ZdGigi.M2.sum=NA,  
                           ZdMomo.M1.sum=NA,ZdMomo.M2.sum=NA,TdKS.M1.sum=NA,TdKS.M2.sum=NA,ZdGigi_4to1.M1.sum=NA,ZdGigi_4to1.M2.sum=NA,  
                           ZdMomo_4to1.M1.sum=NA,ZdMomo_4to1.M2.sum=NA,ZhRIMHU001.M1.sum=NA,ZhRIMHU001.M2.sum=NA,ZmB73.M1.sum=NA,
                           ZmB73.M2.sum=NA,ZmB97.M1.sum=NA,ZmB97.M2.sum=NA,ZmCML103.M1.sum=NA,ZmCML103.M2.sum=NA,
                           ZmCML228.M1.sum=NA,ZmCML228.M2.sum=NA,ZmCML247.M1.sum=NA,ZmCML247.M2.sum=NA,ZmCML277.M1.sum=NA,
                           ZmCML277.M2.sum=NA,ZmCML322.M1.sum=NA,ZmCML322.M2.sum=NA,ZmCML333.M1.sum=NA,ZmCML333.M2.sum=NA,
                           ZmCML52.M1.sum=NA,ZmCML52.M2.sum=NA,ZmCML69.M1.sum=NA,ZmCML69.M2.sum=NA,ZmHP301.M1.sum=NA,
                           ZmHP301.M2.sum=NA,ZmIL14H.M1.sum=NA,ZmIL14H.M2.sum=NA,ZmKi11.M1.sum=NA,ZmKi11.M2.sum=NA,
                           ZmKi3.M1.sum=NA,ZmKi3.M2.sum=NA,ZmKy21.M1.sum=NA,ZmKy21.M2.sum=NA,ZmM162W.M1.sum=NA,
                           ZmM162W.M2.sum=NA,ZmM37W.M1.sum=NA,ZmM37W.M2.sum=NA,ZmMo18W.M1.sum=NA,ZmMo18W.M2.sum=NA,
                           ZmMS71.M1.sum=NA,ZmMS71.M2.sum=NA,ZmNC350.M1.sum=NA,ZmNC350.M2.sum=NA,ZmNC358.M1.sum=NA,
                           ZmNC358.M2.sum=NA,ZmOh43.M1.sum=NA,ZmOh43.M2.sum=NA,ZmOh7b.M1.sum=NA,ZmOh7b.M2.sum=NA,
                           ZmP39.M1.sum=NA,ZmP39.M2.sum=NA,ZmTx303.M1.sum=NA,ZmTx303.M2.sum=NA,ZmTzi8.M1.sum=NA,
                           ZmTzi8.M2.sum=NA,ZnPI615697_4to1.M1.sum=NA,ZnPI615697_4to1.M2.sum=NA,ZnPI615697.M1.sum=NA,ZnPI615697.M2.sum=NA,ZvTIL01.M1.sum=NA,ZvTIL01.M2.sum=NA, 
                           ZvTIL11.M1.sum=NA,ZvTIL11.M2.sum=NA,ZxTIL18.M1.sum=NA,ZxTIL18.M2.sum=NA,ZxTIL25.M1.sum=NA,
                           ZxTIL25.M2.sum=NA)

for(i in 1:length(gene_id)){
  df<-filter(full.fractionation.status, Gene_ID == gene_id[i]) #filter full calls to those for the gene ID
  gene_fractionation[i,"Gene_ID"]<-gene_id[i]
  for(j in 1:length(tripsacinae_genome_IDs)){ #loop through each genome
    #for M1
    if(is.na(all(df[,paste0(tripsacinae_genome_IDs[j],".M1")]))){ #if all the exons for a gene are called NA
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[j],".M1.sum")] <- NA #then gene status is NA
    }else{ #if there is a fractionation call for at least one exon in the gene
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[j],".M1.sum")] <- sum(df[,paste0(tripsacinae_genome_IDs[j],".M1")],na.rm = T)
      #Then sum the calls with na.rm=T
    }
    #for M2
    if(is.na(all(df[,paste0(tripsacinae_genome_IDs[j],".M2")]))){ #if all the exons for a gene are called NA
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[j],".M2.sum")] <- NA #then gene status is NA
    }else{ #if there is a fractionation call for at least one exon in the gene
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[j],".M2.sum")] <- sum(df[,paste0(tripsacinae_genome_IDs[j],".M2")],na.rm = T)
      #Then sum the calls with na.rm=T
    }
  }
}
#There will be a lot of warnings because it has to coerce the numeric to logical to evaluate

#Make a new column to be the status for a given homoeolog
##gene_fractionation[,paste0("TdFL",".M1")]<-NA
#if the sum is greater than 0 and not NA, then give it a 1 (gene is fractionated)
#if the sum is 0 and not NA, then give it a 0 (gene is retained)
#otherwise, it will remain NA
##for(i in 1:nrow(gene_fractionation)){
##  if(gene_fractionation[i,paste0("TdFL",".M1.sum")] > 0 & !is.na(gene_fractionation[i,paste0("TdFL",".M1.sum")])){
##    gene_fractionation[i,paste0("TdFL",".M1")]<-1
##  }else{if(gene_fractionation[i,paste0("TdFL",".M1.sum")] == 0 & !is.na(gene_fractionation[i,paste0("TdFL",".M1.sum")])){
##    gene_fractionation[i,paste0("TdFL",".M1")]<-0
##  }}
##}

for(g in 1:length(tripsacinae_genome_IDs)){
  gene_fractionation[,paste0(tripsacinae_genome_IDs[g],".M1")]<-NA
  for(i in 1:nrow(gene_fractionation)){
    if(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M1.sum")] > 0 & !is.na(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M1.sum")])){
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M1")]<-1 
    }else{if(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M1.sum")] == 0 & !is.na(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M1.sum")])){
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M1")]<-0 
    }}
  }
}

for(g in 1:length(tripsacinae_genome_IDs)){
  gene_fractionation[,paste0(tripsacinae_genome_IDs[g],".M2")]<-NA
  for(i in 1:nrow(gene_fractionation)){
    if(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M2.sum")] > 0 & !is.na(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M2.sum")])){
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M2")]<-1 
    }else{if(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M2.sum")] == 0 & !is.na(gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M2.sum")])){
      gene_fractionation[i,paste0(tripsacinae_genome_IDs[g],".M2")]<-0 
    }}
  }
}

gene_fractionation<-gene_fractionation %>% mutate(TdFL.status = case_when(TdFL.M1 == 0 & TdFL.M2 == 0 ~ "Both_Retained",
                                                                          TdFL.M1 == 1 & TdFL.M2 == 0 ~ "M2_Retained",
                                                                          TdFL.M1 == 0 & TdFL.M2 == 1 ~ "M1_Retained",
                                                                          TdFL.M1 == 1 & TdFL.M2 == 1 ~ "Both_Lost",
                                                                          TdFL.M1 == 1 & is.na(TdFL.M2) ~ "M1_Lost:M2_NA",
                                                                          TdFL.M1 == 0 & is.na(TdFL.M2) ~ "M1_Retained:M2_NA",
                                                                          is.na(TdFL.M1) & TdFL.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                          is.na(TdFL.M1) & TdFL.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                          is.na(TdFL.M1) & is.na(TdFL.M2) ~ "Both_NA"),
                                                  
                                                  TdKS.status = case_when(TdKS.M1 == 0 & TdKS.M2 == 0 ~ "Both_Retained",
                                                                          TdKS.M1 == 1 & TdKS.M2 == 0 ~ "M2_Retained",
                                                                          TdKS.M1 == 0 & TdKS.M2 == 1 ~ "M1_Retained",
                                                                          TdKS.M1 == 1 & TdKS.M2 == 1 ~ "Both_Lost",
                                                                          TdKS.M1 == 1 & is.na(TdKS.M2) ~ "M1_Lost:M2_NA",
                                                                          TdKS.M1 == 0 & is.na(TdKS.M2) ~ "M1_Retained:M2_NA",
                                                                          is.na(TdKS.M1) & TdKS.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                          is.na(TdKS.M1) & TdKS.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                          is.na(TdKS.M1) & is.na(TdKS.M2) ~ "Both_NA"),
                                                  
                                                  ZdGigi_4to1.status = case_when(ZdGigi_4to1.M1 == 0 & ZdGigi_4to1.M2 == 0 ~ "Both_Retained",
                                                                                 ZdGigi_4to1.M1 == 1 & ZdGigi_4to1.M2 == 0 ~ "M2_Retained",
                                                                          ZdGigi_4to1.M1 == 0 & ZdGigi_4to1.M2 == 1 ~ "M1_Retained",
                                                                          ZdGigi_4to1.M1 == 1 & ZdGigi_4to1.M2 == 1 ~ "Both_Lost",
                                                                          ZdGigi_4to1.M1 == 1 & is.na(ZdGigi_4to1.M2) ~ "M1_Lost:M2_NA",
                                                                          ZdGigi_4to1.M1 == 0 & is.na(ZdGigi_4to1.M2) ~ "M1_Retained:M2_NA",
                                                                          is.na(ZdGigi_4to1.M1) & ZdGigi_4to1.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                          is.na(ZdGigi_4to1.M1) & ZdGigi_4to1.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                          is.na(ZdGigi_4to1.M1) & is.na(ZdGigi_4to1.M2) ~ "Both_NA"),
                                                  
                                                  ZdMomo_4to1.status = case_when(ZdMomo_4to1.M1 == 0 & ZdMomo_4to1.M2 == 0 ~ "Both_Retained",
                                                                                 ZdMomo_4to1.M1 == 1 & ZdMomo_4to1.M2 == 0 ~ "M2_Retained",
                                                                                 ZdMomo_4to1.M1 == 0 & ZdMomo_4to1.M2 == 1 ~ "M1_Retained",
                                                                                 ZdMomo_4to1.M1 == 1 & ZdMomo_4to1.M2 == 1 ~ "Both_Lost",
                                                                                 ZdMomo_4to1.M1 == 1 & is.na(ZdMomo_4to1.M2) ~ "M1_Lost:M2_NA",
                                                                                 ZdMomo_4to1.M1 == 0 & is.na(ZdMomo_4to1.M2) ~ "M1_Retained:M2_NA",
                                                                          is.na(ZdMomo_4to1.M1) & ZdMomo_4to1.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                          is.na(ZdMomo_4to1.M1) & ZdMomo_4to1.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                          is.na(ZdMomo_4to1.M1) & is.na(ZdMomo_4to1.M2) ~ "Both_NA"),
                                                  
                                                  ZnPI615697_4to1.status = case_when(ZnPI615697_4to1.M1 == 0 & ZnPI615697_4to1.M2 == 0 ~ "Both_Retained",
                                                                                     ZnPI615697_4to1.M1 == 1 & ZnPI615697_4to1.M2 == 0 ~ "M2_Retained",
                                                                                     ZnPI615697_4to1.M1 == 0 & ZnPI615697_4to1.M2 == 1 ~ "M1_Retained",
                                                                                     ZnPI615697_4to1.M1 == 1 & ZnPI615697_4to1.M2 == 1 ~ "Both_Lost",
                                                                                     ZnPI615697_4to1.M1 == 1 & is.na(ZnPI615697_4to1.M2) ~ "M1_Lost:M2_NA",
                                                                                     ZnPI615697_4to1.M1 == 0 & is.na(ZnPI615697_4to1.M2) ~ "M1_Retained:M2_NA",
                                                                          is.na(ZnPI615697_4to1.M1) & ZnPI615697_4to1.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                          is.na(ZnPI615697_4to1.M1) & ZnPI615697_4to1.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                          is.na(ZnPI615697_4to1.M1) & is.na(ZnPI615697_4to1.M2) ~ "Both_NA"),
                                                  
                                                  ZdGigi.status = case_when(ZdGigi.M1 == 0 & ZdGigi.M2 == 0 ~ "Both_Retained",
                                                                            ZdGigi.M1 == 1 & ZdGigi.M2 == 0 ~ "M2_Retained",
                                                                            ZdGigi.M1 == 0 & ZdGigi.M2 == 1 ~ "M1_Retained",
                                                                            ZdGigi.M1 == 1 & ZdGigi.M2 == 1 ~ "Both_Lost",
                                                                            ZdGigi.M1 == 1 & is.na(ZdGigi.M2) ~ "M1_Lost:M2_NA",
                                                                            ZdGigi.M1 == 0 & is.na(ZdGigi.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZdGigi.M1) & ZdGigi.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZdGigi.M1) & ZdGigi.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZdGigi.M1) & is.na(ZdGigi.M2) ~ "Both_NA"),
                                                  
                                                  ZdMomo.status = case_when(ZdMomo.M1 == 0 & ZdMomo.M2 == 0 ~ "Both_Retained",
                                                                            ZdMomo.M1 == 1 & ZdMomo.M2 == 0 ~ "M2_Retained",
                                                                            ZdMomo.M1 == 0 & ZdMomo.M2 == 1 ~ "M1_Retained",
                                                                            ZdMomo.M1 == 1 & ZdMomo.M2 == 1 ~ "Both_Lost",
                                                                            ZdMomo.M1 == 1 & is.na(ZdMomo.M2) ~ "M1_Lost:M2_NA",
                                                                            ZdMomo.M1 == 0 & is.na(ZdMomo.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZdMomo.M1) & ZdMomo.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZdMomo.M1) & ZdMomo.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZdMomo.M1) & is.na(ZdMomo.M2) ~ "Both_NA"),
                                                  
                                                  ZhRIMHU001.status = case_when(ZhRIMHU001.M1 == 0 & ZhRIMHU001.M2 == 0 ~ "Both_Retained",
                                                                                ZhRIMHU001.M1 == 1 & ZhRIMHU001.M2 == 0 ~ "M2_Retained",
                                                                                ZhRIMHU001.M1 == 0 & ZhRIMHU001.M2 == 1 ~ "M1_Retained",
                                                                                ZhRIMHU001.M1 == 1 & ZhRIMHU001.M2 == 1 ~ "Both_Lost",
                                                                                ZhRIMHU001.M1 == 1 & is.na(ZhRIMHU001.M2) ~ "M1_Lost:M2_NA",
                                                                                ZhRIMHU001.M1 == 0 & is.na(ZhRIMHU001.M2) ~ "M1_Retained:M2_NA",
                                                                                is.na(ZhRIMHU001.M1) & ZhRIMHU001.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                is.na(ZhRIMHU001.M1) & ZhRIMHU001.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                is.na(ZhRIMHU001.M1) & is.na(ZhRIMHU001.M2) ~ "Both_NA"),
                                                  
                                                  ZmB73.status = case_when(ZmB73.M1 == 0 & ZmB73.M2 == 0 ~ "Both_Retained",
                                                                           ZmB73.M1 == 1 & ZmB73.M2 == 0 ~ "M2_Retained",
                                                                           ZmB73.M1 == 0 & ZmB73.M2 == 1 ~ "M1_Retained",
                                                                           ZmB73.M1 == 1 & ZmB73.M2 == 1 ~ "Both_Lost",
                                                                           ZmB73.M1 == 1 & is.na(ZmB73.M2) ~ "M1_Lost:M2_NA",
                                                                           ZmB73.M1 == 0 & is.na(ZmB73.M2) ~ "M1_Retained:M2_NA",
                                                                           is.na(ZmB73.M1) & ZmB73.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                           is.na(ZmB73.M1) & ZmB73.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                           is.na(ZmB73.M1) & is.na(ZmB73.M2) ~ "Both_NA"),
                                                  
                                                  ZmB97.status = case_when(ZmB97.M1 == 0 & ZmB97.M2 == 0 ~ "Both_Retained",
                                                                           ZmB97.M1 == 1 & ZmB97.M2 == 0 ~ "M2_Retained",
                                                                           ZmB97.M1 == 0 & ZmB97.M2 == 1 ~ "M1_Retained",
                                                                           ZmB97.M1 == 1 & ZmB97.M2 == 1 ~ "Both_Lost",
                                                                           ZmB97.M1 == 1 & is.na(ZmB97.M2) ~ "M1_Lost:M2_NA",
                                                                           ZmB97.M1 == 0 & is.na(ZmB97.M2) ~ "M1_Retained:M2_NA",
                                                                           is.na(ZmB97.M1) & ZmB97.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                           is.na(ZmB97.M1) & ZmB97.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                           is.na(ZmB97.M1) & is.na(ZmB97.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML103.status = case_when(ZmCML103.M1 == 0 & ZmCML103.M2 == 0 ~ "Both_Retained",
                                                                              ZmCML103.M1 == 1 & ZmCML103.M2 == 0 ~ "M2_Retained",
                                                                              ZmCML103.M1 == 0 & ZmCML103.M2 == 1 ~ "M1_Retained",
                                                                              ZmCML103.M1 == 1 & ZmCML103.M2 == 1 ~ "Both_Lost",
                                                                              ZmCML103.M1 == 1 & is.na(ZmCML103.M2) ~ "M1_Lost:M2_NA",
                                                                              ZmCML103.M1 == 0 & is.na(ZmCML103.M2) ~ "M1_Retained:M2_NA",
                                                                              is.na(ZmCML103.M1) & ZmCML103.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                              is.na(ZmCML103.M1) & ZmCML103.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                              is.na(ZmCML103.M1) & is.na(ZmCML103.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML228.status = case_when(ZmCML228.M1 == 0 & ZmCML228.M2 == 0 ~ "Both_Retained",
                                                                              ZmCML228.M1 == 1 & ZmCML228.M2 == 0 ~ "M2_Retained",
                                                                              ZmCML228.M1 == 0 & ZmCML228.M2 == 1 ~ "M1_Retained",
                                                                              ZmCML228.M1 == 1 & ZmCML228.M2 == 1 ~ "Both_Lost",
                                                                              ZmCML228.M1 == 1 & is.na(ZmCML228.M2) ~ "M1_Lost:M2_NA",
                                                                              ZmCML228.M1 == 0 & is.na(ZmCML228.M2) ~ "M1_Retained:M2_NA",
                                                                              is.na(ZmCML228.M1) & ZmCML228.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                              is.na(ZmCML228.M1) & ZmCML228.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                              is.na(ZmCML228.M1) & is.na(ZmCML228.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML247.status = case_when(ZmCML247.M1 == 0 & ZmCML247.M2 == 0 ~ "Both_Retained",
                                                                              ZmCML247.M1 == 1 & ZmCML247.M2 == 0 ~ "M2_Retained",
                                                                              ZmCML247.M1 == 0 & ZmCML247.M2 == 1 ~ "M1_Retained",
                                                                              ZmCML247.M1 == 1 & ZmCML247.M2 == 1 ~ "Both_Lost",
                                                                              ZmCML247.M1 == 1 & is.na(ZmCML247.M2) ~ "M1_Lost:M2_NA",
                                                                              ZmCML247.M1 == 0 & is.na(ZmCML247.M2) ~ "M1_Retained:M2_NA",
                                                                              is.na(ZmCML247.M1) & ZmCML247.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                              is.na(ZmCML247.M1) & ZmCML247.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                              is.na(ZmCML247.M1) & is.na(ZmCML247.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML277.status = case_when(ZmCML277.M1 == 0 & ZmCML277.M2 == 0 ~ "Both_Retained",
                                                                              ZmCML277.M1 == 1 & ZmCML277.M2 == 0 ~ "M2_Retained",
                                                                              ZmCML277.M1 == 0 & ZmCML277.M2 == 1 ~ "M1_Retained",
                                                                              ZmCML277.M1 == 1 & ZmCML277.M2 == 1 ~ "Both_Lost",
                                                                              ZmCML277.M1 == 1 & is.na(ZmCML277.M2) ~ "M1_Lost:M2_NA",
                                                                              ZmCML277.M1 == 0 & is.na(ZmCML277.M2) ~ "M1_Retained:M2_NA",
                                                                              is.na(ZmCML277.M1) & ZmCML277.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                              is.na(ZmCML277.M1) & ZmCML277.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                              is.na(ZmCML277.M1) & is.na(ZmCML277.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML322.status = case_when(ZmCML322.M1 == 0 & ZmCML322.M2 == 0 ~ "Both_Retained",
                                                                              ZmCML322.M1 == 1 & ZmCML322.M2 == 0 ~ "M2_Retained",
                                                                              ZmCML322.M1 == 0 & ZmCML322.M2 == 1 ~ "M1_Retained",
                                                                              ZmCML322.M1 == 1 & ZmCML322.M2 == 1 ~ "Both_Lost",
                                                                              ZmCML322.M1 == 1 & is.na(ZmCML322.M2) ~ "M1_Lost:M2_NA",
                                                                              ZmCML322.M1 == 0 & is.na(ZmCML322.M2) ~ "M1_Retained:M2_NA",
                                                                              is.na(ZmCML322.M1) & ZmCML322.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                              is.na(ZmCML322.M1) & ZmCML322.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                              is.na(ZmCML322.M1) & is.na(ZmCML322.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML333.status = case_when(ZmCML333.M1 == 0 & ZmCML333.M2 == 0 ~ "Both_Retained",
                                                                              ZmCML333.M1 == 1 & ZmCML333.M2 == 0 ~ "M2_Retained",
                                                                              ZmCML333.M1 == 0 & ZmCML333.M2 == 1 ~ "M1_Retained",
                                                                              ZmCML333.M1 == 1 & ZmCML333.M2 == 1 ~ "Both_Lost",
                                                                              ZmCML333.M1 == 1 & is.na(ZmCML333.M2) ~ "M1_Lost:M2_NA",
                                                                              ZmCML333.M1 == 0 & is.na(ZmCML333.M2) ~ "M1_Retained:M2_NA",
                                                                              is.na(ZmCML333.M1) & ZmCML333.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                              is.na(ZmCML333.M1) & ZmCML333.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                              is.na(ZmCML333.M1) & is.na(ZmCML333.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML52.status = case_when(ZmCML52.M1 == 0 & ZmCML52.M2 == 0 ~ "Both_Retained",
                                                                             ZmCML52.M1 == 1 & ZmCML52.M2 == 0 ~ "M2_Retained",
                                                                             ZmCML52.M1 == 0 & ZmCML52.M2 == 1 ~ "M1_Retained",
                                                                             ZmCML52.M1 == 1 & ZmCML52.M2 == 1 ~ "Both_Lost",
                                                                             ZmCML52.M1 == 1 & is.na(ZmCML52.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmCML52.M1 == 0 & is.na(ZmCML52.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmCML52.M1) & ZmCML52.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmCML52.M1) & ZmCML52.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmCML52.M1) & is.na(ZmCML52.M2) ~ "Both_NA"),
                                                  
                                                  ZmCML69.status = case_when(ZmCML69.M1 == 0 & ZmCML69.M2 == 0 ~ "Both_Retained",
                                                                             ZmCML69.M1 == 1 & ZmCML69.M2 == 0 ~ "M2_Retained",
                                                                             ZmCML69.M1 == 0 & ZmCML69.M2 == 1 ~ "M1_Retained",
                                                                             ZmCML69.M1 == 1 & ZmCML69.M2 == 1 ~ "Both_Lost",
                                                                             ZmCML69.M1 == 1 & is.na(ZmCML69.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmCML69.M1 == 0 & is.na(ZmCML69.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmCML69.M1) & ZmCML69.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmCML69.M1) & ZmCML69.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmCML69.M1) & is.na(ZmCML69.M2) ~ "Both_NA"),
                                                  
                                                  ZmHP301.status = case_when(ZmHP301.M1 == 0 & ZmHP301.M2 == 0 ~ "Both_Retained",
                                                                             ZmHP301.M1 == 1 & ZmHP301.M2 == 0 ~ "M2_Retained",
                                                                             ZmHP301.M1 == 0 & ZmHP301.M2 == 1 ~ "M1_Retained",
                                                                             ZmHP301.M1 == 1 & ZmHP301.M2 == 1 ~ "Both_Lost",
                                                                             ZmHP301.M1 == 1 & is.na(ZmHP301.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmHP301.M1 == 0 & is.na(ZmHP301.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmHP301.M1) & ZmHP301.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmHP301.M1) & ZmHP301.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmHP301.M1) & is.na(ZmHP301.M2) ~ "Both_NA"),
                                                  
                                                  ZmIL14H.status = case_when(ZmIL14H.M1 == 0 & ZmIL14H.M2 == 0 ~ "Both_Retained",
                                                                             ZmIL14H.M1 == 1 & ZmIL14H.M2 == 0 ~ "M2_Retained",
                                                                             ZmIL14H.M1 == 0 & ZmIL14H.M2 == 1 ~ "M1_Retained",
                                                                             ZmIL14H.M1 == 1 & ZmIL14H.M2 == 1 ~ "Both_Lost",
                                                                             ZmIL14H.M1 == 1 & is.na(ZmIL14H.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmIL14H.M1 == 0 & is.na(ZmIL14H.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmIL14H.M1) & ZmIL14H.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmIL14H.M1) & ZmIL14H.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmIL14H.M1) & is.na(ZmIL14H.M2) ~ "Both_NA"),
                                                  
                                                  ZmKi11.status = case_when(ZmKi11.M1 == 0 & ZmKi11.M2 == 0 ~ "Both_Retained",
                                                                            ZmKi11.M1 == 1 & ZmKi11.M2 == 0 ~ "M2_Retained",
                                                                            ZmKi11.M1 == 0 & ZmKi11.M2 == 1 ~ "M1_Retained",
                                                                            ZmKi11.M1 == 1 & ZmKi11.M2 == 1 ~ "Both_Lost",
                                                                            ZmKi11.M1 == 1 & is.na(ZmKi11.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmKi11.M1 == 0 & is.na(ZmKi11.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmKi11.M1) & ZmKi11.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmKi11.M1) & ZmKi11.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmKi11.M1) & is.na(ZmKi11.M2) ~ "Both_NA"),
                                                  
                                                  ZmKi3.status = case_when(ZmKi3.M1 == 0 & ZmKi3.M2 == 0 ~ "Both_Retained",
                                                                           ZmKi3.M1 == 1 & ZmKi3.M2 == 0 ~ "M2_Retained",
                                                                           ZmKi3.M1 == 0 & ZmKi3.M2 == 1 ~ "M1_Retained",
                                                                           ZmKi3.M1 == 1 & ZmKi3.M2 == 1 ~ "Both_Lost",
                                                                           ZmKi3.M1 == 1 & is.na(ZmKi3.M2) ~ "M1_Lost:M2_NA",
                                                                           ZmKi3.M1 == 0 & is.na(ZmKi3.M2) ~ "M1_Retained:M2_NA",
                                                                           is.na(ZmKi3.M1) & ZmKi3.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                           is.na(ZmKi3.M1) & ZmKi3.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                           is.na(ZmKi3.M1) & is.na(ZmKi3.M2) ~ "Both_NA"),
                                                  
                                                  ZmKy21.status = case_when(ZmKy21.M1 == 0 & ZmKy21.M2 == 0 ~ "Both_Retained",
                                                                            ZmKy21.M1 == 1 & ZmKy21.M2 == 0 ~ "M2_Retained",
                                                                            ZmKy21.M1 == 0 & ZmKy21.M2 == 1 ~ "M1_Retained",
                                                                            ZmKy21.M1 == 1 & ZmKy21.M2 == 1 ~ "Both_Lost",
                                                                            ZmKy21.M1 == 1 & is.na(ZmKy21.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmKy21.M1 == 0 & is.na(ZmKy21.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmKy21.M1) & ZmKy21.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmKy21.M1) & ZmKy21.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmKy21.M1) & is.na(ZmKy21.M2) ~ "Both_NA"),
                                                  
                                                  ZmM162W.status = case_when(ZmM162W.M1 == 0 & ZmM162W.M2 == 0 ~ "Both_Retained",
                                                                             ZmM162W.M1 == 1 & ZmM162W.M2 == 0 ~ "M2_Retained",
                                                                             ZmM162W.M1 == 0 & ZmM162W.M2 == 1 ~ "M1_Retained",
                                                                             ZmM162W.M1 == 1 & ZmM162W.M2 == 1 ~ "Both_Lost",
                                                                             ZmM162W.M1 == 1 & is.na(ZmM162W.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmM162W.M1 == 0 & is.na(ZmM162W.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmM162W.M1) & ZmM162W.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmM162W.M1) & ZmM162W.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmM162W.M1) & is.na(ZmM162W.M2) ~ "Both_NA"),
                                                  
                                                  ZmM37W.status = case_when(ZmM37W.M1 == 0 & ZmM37W.M2 == 0 ~ "Both_Retained",
                                                                            ZmM37W.M1 == 1 & ZmM37W.M2 == 0 ~ "M2_Retained",
                                                                            ZmM37W.M1 == 0 & ZmM37W.M2 == 1 ~ "M1_Retained",
                                                                            ZmM37W.M1 == 1 & ZmM37W.M2 == 1 ~ "Both_Lost",
                                                                            ZmM37W.M1 == 1 & is.na(ZmM37W.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmM37W.M1 == 0 & is.na(ZmM37W.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmM37W.M1) & ZmM37W.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmM37W.M1) & ZmM37W.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmM37W.M1) & is.na(ZmM37W.M2) ~ "Both_NA"),
                                                  
                                                  ZmMo18W.status = case_when(ZmMo18W.M1 == 0 & ZmMo18W.M2 == 0 ~ "Both_Retained",
                                                                             ZmMo18W.M1 == 1 & ZmMo18W.M2 == 0 ~ "M2_Retained",
                                                                             ZmMo18W.M1 == 0 & ZmMo18W.M2 == 1 ~ "M1_Retained",
                                                                             ZmMo18W.M1 == 1 & ZmMo18W.M2 == 1 ~ "Both_Lost",
                                                                             ZmMo18W.M1 == 1 & is.na(ZmMo18W.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmMo18W.M1 == 0 & is.na(ZmMo18W.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmMo18W.M1) & ZmMo18W.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmMo18W.M1) & ZmMo18W.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmMo18W.M1) & is.na(ZmMo18W.M2) ~ "Both_NA"),
                                                  
                                                  ZmMS71.status = case_when(ZmMS71.M1 == 0 & ZmMS71.M2 == 0 ~ "Both_Retained",
                                                                            ZmMS71.M1 == 1 & ZmMS71.M2 == 0 ~ "M2_Retained",
                                                                            ZmMS71.M1 == 0 & ZmMS71.M2 == 1 ~ "M1_Retained",
                                                                            ZmMS71.M1 == 1 & ZmMS71.M2 == 1 ~ "Both_Lost",
                                                                            ZmMS71.M1 == 1 & is.na(ZmMS71.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmMS71.M1 == 0 & is.na(ZmMS71.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmMS71.M1) & ZmMS71.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmMS71.M1) & ZmMS71.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmMS71.M1) & is.na(ZmMS71.M2) ~ "Both_NA"),
                                                  
                                                  ZmNC350.status = case_when(ZmNC350.M1 == 0 & ZmNC350.M2 == 0 ~ "Both_Retained",
                                                                             ZmNC350.M1 == 1 & ZmNC350.M2 == 0 ~ "M2_Retained",
                                                                             ZmNC350.M1 == 0 & ZmNC350.M2 == 1 ~ "M1_Retained",
                                                                             ZmNC350.M1 == 1 & ZmNC350.M2 == 1 ~ "Both_Lost",
                                                                             ZmNC350.M1 == 1 & is.na(ZmNC350.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmNC350.M1 == 0 & is.na(ZmNC350.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmNC350.M1) & ZmNC350.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmNC350.M1) & ZmNC350.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmNC350.M1) & is.na(ZmNC350.M2) ~ "Both_NA"),
                                                  
                                                  ZmNC358.status = case_when(ZmNC358.M1 == 0 & ZmNC358.M2 == 0 ~ "Both_Retained",
                                                                             ZmNC358.M1 == 1 & ZmNC358.M2 == 0 ~ "M2_Retained",
                                                                             ZmNC358.M1 == 0 & ZmNC358.M2 == 1 ~ "M1_Retained",
                                                                             ZmNC358.M1 == 1 & ZmNC358.M2 == 1 ~ "Both_Lost",
                                                                             ZmNC358.M1 == 1 & is.na(ZmNC358.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmNC358.M1 == 0 & is.na(ZmNC358.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmNC358.M1) & ZmNC358.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmNC358.M1) & ZmNC358.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmNC358.M1) & is.na(ZmNC358.M2) ~ "Both_NA"),
                                                  
                                                  ZmOh43.status = case_when(ZmOh43.M1 == 0 & ZmOh43.M2 == 0 ~ "Both_Retained",
                                                                            ZmOh43.M1 == 1 & ZmOh43.M2 == 0 ~ "M2_Retained",
                                                                            ZmOh43.M1 == 0 & ZmOh43.M2 == 1 ~ "M1_Retained",
                                                                            ZmOh43.M1 == 1 & ZmOh43.M2 == 1 ~ "Both_Lost",
                                                                            ZmOh43.M1 == 1 & is.na(ZmOh43.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmOh43.M1 == 0 & is.na(ZmOh43.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmOh43.M1) & ZmOh43.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmOh43.M1) & ZmOh43.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmOh43.M1) & is.na(ZmOh43.M2) ~ "Both_NA"),
                                                  
                                                  ZmOh7b.status = case_when(ZmOh7b.M1 == 0 & ZmOh7b.M2 == 0 ~ "Both_Retained",
                                                                            ZmOh7b.M1 == 1 & ZmOh7b.M2 == 0 ~ "M2_Retained",
                                                                            ZmOh7b.M1 == 0 & ZmOh7b.M2 == 1 ~ "M1_Retained",
                                                                            ZmOh7b.M1 == 1 & ZmOh7b.M2 == 1 ~ "Both_Lost",
                                                                            ZmOh7b.M1 == 1 & is.na(ZmOh7b.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmOh7b.M1 == 0 & is.na(ZmOh7b.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmOh7b.M1) & ZmOh7b.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmOh7b.M1) & ZmOh7b.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmOh7b.M1) & is.na(ZmOh7b.M2) ~ "Both_NA"),
                                                  
                                                  ZmP39.status = case_when(ZmP39.M1 == 0 & ZmP39.M2 == 0 ~ "Both_Retained",
                                                                           ZmP39.M1 == 1 & ZmP39.M2 == 0 ~ "M2_Retained",
                                                                           ZmP39.M1 == 0 & ZmP39.M2 == 1 ~ "M1_Retained",
                                                                           ZmP39.M1 == 1 & ZmP39.M2 == 1 ~ "Both_Lost",
                                                                           ZmP39.M1 == 1 & is.na(ZmP39.M2) ~ "M1_Lost:M2_NA",
                                                                           ZmP39.M1 == 0 & is.na(ZmP39.M2) ~ "M1_Retained:M2_NA",
                                                                           is.na(ZmP39.M1) & ZmP39.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                           is.na(ZmP39.M1) & ZmP39.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                           is.na(ZmP39.M1) & is.na(ZmP39.M2) ~ "Both_NA"),
                                                  
                                                  ZmTx303.status = case_when(ZmTx303.M1 == 0 & ZmTx303.M2 == 0 ~ "Both_Retained",
                                                                             ZmTx303.M1 == 1 & ZmTx303.M2 == 0 ~ "M2_Retained",
                                                                             ZmTx303.M1 == 0 & ZmTx303.M2 == 1 ~ "M1_Retained",
                                                                             ZmTx303.M1 == 1 & ZmTx303.M2 == 1 ~ "Both_Lost",
                                                                             ZmTx303.M1 == 1 & is.na(ZmTx303.M2) ~ "M1_Lost:M2_NA",
                                                                             ZmTx303.M1 == 0 & is.na(ZmTx303.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZmTx303.M1) & ZmTx303.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZmTx303.M1) & ZmTx303.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZmTx303.M1) & is.na(ZmTx303.M2) ~ "Both_NA"),
                                                  
                                                  ZmTzi8.status = case_when(ZmTzi8.M1 == 0 & ZmTzi8.M2 == 0 ~ "Both_Retained",
                                                                            ZmTzi8.M1 == 1 & ZmTzi8.M2 == 0 ~ "M2_Retained",
                                                                            ZmTzi8.M1 == 0 & ZmTzi8.M2 == 1 ~ "M1_Retained",
                                                                            ZmTzi8.M1 == 1 & ZmTzi8.M2 == 1 ~ "Both_Lost",
                                                                            ZmTzi8.M1 == 1 & is.na(ZmTzi8.M2) ~ "M1_Lost:M2_NA",
                                                                            ZmTzi8.M1 == 0 & is.na(ZmTzi8.M2) ~ "M1_Retained:M2_NA",
                                                                            is.na(ZmTzi8.M1) & ZmTzi8.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                            is.na(ZmTzi8.M1) & ZmTzi8.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                            is.na(ZmTzi8.M1) & is.na(ZmTzi8.M2) ~ "Both_NA"),
                                                  
                                                  ZnPI615697.status = case_when(ZnPI615697.M1 == 0 & ZnPI615697.M2 == 0 ~ "Both_Retained",
                                                                                ZnPI615697.M1 == 1 & ZnPI615697.M2 == 0 ~ "M2_Retained",
                                                                                ZnPI615697.M1 == 0 & ZnPI615697.M2 == 1 ~ "M1_Retained",
                                                                                ZnPI615697.M1 == 1 & ZnPI615697.M2 == 1 ~ "Both_Lost",
                                                                                ZnPI615697.M1 == 1 & is.na(ZnPI615697.M2) ~ "M1_Lost:M2_NA",
                                                                                ZnPI615697.M1 == 0 & is.na(ZnPI615697.M2) ~ "M1_Retained:M2_NA",
                                                                                is.na(ZnPI615697.M1) & ZnPI615697.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                                is.na(ZnPI615697.M1) & ZnPI615697.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                                is.na(ZnPI615697.M1) & is.na(ZnPI615697.M2) ~ "Both_NA"),
                                                  
                                                  ZvTIL01.status = case_when(ZvTIL01.M1 == 0 & ZvTIL01.M2 == 0 ~ "Both_Retained",
                                                                             ZvTIL01.M1 == 1 & ZvTIL01.M2 == 0 ~ "M2_Retained",
                                                                             ZvTIL01.M1 == 0 & ZvTIL01.M2 == 1 ~ "M1_Retained",
                                                                             ZvTIL01.M1 == 1 & ZvTIL01.M2 == 1 ~ "Both_Lost",
                                                                             ZvTIL01.M1 == 1 & is.na(ZvTIL01.M2) ~ "M1_Lost:M2_NA",
                                                                             ZvTIL01.M1 == 0 & is.na(ZvTIL01.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZvTIL01.M1) & ZvTIL01.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZvTIL01.M1) & ZvTIL01.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZvTIL01.M1) & is.na(ZvTIL01.M2) ~ "Both_NA"),
                                                  
                                                  ZvTIL11.status = case_when(ZvTIL11.M1 == 0 & ZvTIL11.M2 == 0 ~ "Both_Retained",
                                                                             ZvTIL11.M1 == 1 & ZvTIL11.M2 == 0 ~ "M2_Retained",
                                                                             ZvTIL11.M1 == 0 & ZvTIL11.M2 == 1 ~ "M1_Retained",
                                                                             ZvTIL11.M1 == 1 & ZvTIL11.M2 == 1 ~ "Both_Lost",
                                                                             ZvTIL11.M1 == 1 & is.na(ZvTIL11.M2) ~ "M1_Lost:M2_NA",
                                                                             ZvTIL11.M1 == 0 & is.na(ZvTIL11.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZvTIL11.M1) & ZvTIL11.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZvTIL11.M1) & ZvTIL11.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZvTIL11.M1) & is.na(ZvTIL11.M2) ~ "Both_NA"),
                                                  
                                                  ZxTIL18.status = case_when(ZxTIL18.M1 == 0 & ZxTIL18.M2 == 0 ~ "Both_Retained",
                                                                             ZxTIL18.M1 == 1 & ZxTIL18.M2 == 0 ~ "M2_Retained",
                                                                             ZxTIL18.M1 == 0 & ZxTIL18.M2 == 1 ~ "M1_Retained",
                                                                             ZxTIL18.M1 == 1 & ZxTIL18.M2 == 1 ~ "Both_Lost",
                                                                             ZxTIL18.M1 == 1 & is.na(ZxTIL18.M2) ~ "M1_Lost:M2_NA",
                                                                             ZxTIL18.M1 == 0 & is.na(ZxTIL18.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZxTIL18.M1) & ZxTIL18.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZxTIL18.M1) & ZxTIL18.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZxTIL18.M1) & is.na(ZxTIL18.M2) ~ "Both_NA"),
                                                  
                                                  ZxTIL25.status = case_when(ZxTIL25.M1 == 0 & ZxTIL25.M2 == 0 ~ "Both_Retained",
                                                                             ZxTIL25.M1 == 1 & ZxTIL25.M2 == 0 ~ "M2_Retained",
                                                                             ZxTIL25.M1 == 0 & ZxTIL25.M2 == 1 ~ "M1_Retained",
                                                                             ZxTIL25.M1 == 1 & ZxTIL25.M2 == 1 ~ "Both_Lost",
                                                                             ZxTIL25.M1 == 1 & is.na(ZxTIL25.M2) ~ "M1_Lost:M2_NA",
                                                                             ZxTIL25.M1 == 0 & is.na(ZxTIL25.M2) ~ "M1_Retained:M2_NA",
                                                                             is.na(ZxTIL25.M1) & ZxTIL25.M2 == 1 ~ "M1_NA:M2_Lost",
                                                                             is.na(ZxTIL25.M1) & ZxTIL25.M2 == 0 ~ "M1_NA:M2_Retained",
                                                                             is.na(ZxTIL25.M1) & is.na(ZxTIL25.M2) ~ "Both_NA"),
                                                  
)

gene_fractionation %>% select(-c(contains("sum"),ends_with("M1"),ends_with("M2"))) %>%
  pivot_longer(cols = contains(".status"), names_to = "Subgenome", values_to = "Status")%>%
  mutate(Genome = str_split(Subgenome,"[.]",simplify=T)[,1],
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")),
         Status = factor(Status, levels = c("Both_Retained","M1_Retained","M2_Retained", "Both_Lost","M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA"))) %>%
  na.omit()%>%
  ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  xlab("Number of Genes")+ylab("")+
  ggtitle("Fractionation Status by Gene")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByGene.AllGenomes.png", device="png",dpi=500,width = 8, height = 7)

gene_fractionation %>% select(-c(contains("sum"),ends_with("M1"),ends_with("M2"))) %>%
  pivot_longer(cols = contains(".status"), names_to = "Subgenome", values_to = "Status")%>%
  mutate(Genome = str_split(Subgenome,"[.]",simplify=T)[,1],
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","TdFL","TdKS")),
         Status = factor(Status, levels = c("Both_Retained","M1_Retained","M2_Retained", "Both_Lost","M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA"))) %>%
  na.omit()%>% filter(!Genome %in% c("ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1"))%>%
  ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  xlab("Number of Genes")+ylab("")+
  ggtitle("Fractionation Status by Gene")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByGene.noZdZn.png", device="png",dpi=500,width = 8, height = 7)

gene_fractionation%>% select(-c(contains("sum"),ends_with("M1"),ends_with("M2"))) %>%
  pivot_longer(cols = contains(".status"), names_to = "Subgenome", values_to = "Status")%>%
  mutate(Genome = str_split(Subgenome,"[.]",simplify=T)[,1],
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","TdFL"))) %>%
  na.omit()%>%
  filter(Status == "Both_Retained") %>% group_by(Genome) %>% count() %>% as_tibble() %>% 
  reframe(mean=mean(n),min=min(n),max=max(n))
#both homoeologs retained mean=621.   min=304   max=857

gene_fractionation%>% select(-c(contains("sum"),ends_with("M1"),ends_with("M2"))) %>%
  pivot_longer(cols = contains(".status"), names_to = "Subgenome", values_to = "Status")%>%
  mutate(Genome = str_split(Subgenome,"[.]",simplify=T)[,1],
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","TdFL"))) %>%
  na.omit()%>%
  filter(Status %in% c("M1_Retained","M2_Retained")) %>% group_by(Genome) %>% count() %>% as_tibble() %>% 
  reframe(mean=mean(n),min=min(n),max=max(n))
#M1 or M2 homoeolog retained as singleton: mean=2423.  min=1256  max=2556

gene_fractionation%>% select(-c(contains("sum"),ends_with("M1"),ends_with("M2"))) %>%
  pivot_longer(cols = contains(".status"), names_to = "Subgenome", values_to = "Status")%>%
  mutate(Genome = str_split(Subgenome,"[.]",simplify=T)[,1],
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))) %>%
  na.omit()%>%
  filter(Status == "Both_Lost") %>% group_by(Genome) %>% count() %>% as_tibble() %>% reframe(mean=mean(n),min=min(n),max=max(n))
#both homoeologs lost: mean =7478, min=3945, max=7861

gene_fractionation %>% select(c("Gene_ID",contains("status"))) %>% 
  pivot_longer(cols = ends_with("status"), names_to = "Genome", values_to = "Status" ) %>%
  mutate(Genome = str_remove(Genome,".status"),
         Genome = factor(Genome, levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                          "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                          "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                          "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")),
         Status = factor(Status, levels = c("Both_Retained","M1_Retained","M2_Retained", "Both_Lost","M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA"))) %>%
  group_by(Genome,Status) %>% count() %>% 
  mutate(Percent = (n/nrow(gene_fractionation))*100) %>%
  ggplot(aes(y=Status,x=Percent))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  xlab("Percentage of Genes")+ylab("")+
  ggtitle("Fractionation Status by Gene")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByPercentage.AllGenomes.png", device="png",dpi=500,width = 8, height = 7)


#How much differential fractionation is there? 
#Where for a given exon, the status is M1_Retained in one genome, but the in another genome it's M2_Retained
#this object may have been made in the spatial fractionation section
long.full.fractionation.status %>% group_by(ID) %>% count(Status)

#How much missingness is there across subgenomes by exon?
subgenome_NA_Counts<-lapply(select(full.fractionation.status, c(ends_with(".M1"),ends_with(".M2"))), function(x) sum(is.na(x))) %>% 
  as_tibble() %>% pivot_longer(cols = everything(), values_to = "NA_Count") %>% 
  mutate(Subgenome = str_split(name, "[.]", simplify = T)[,2],
         Genome = str_split(name, "[.]", simplify = T)[,1] %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                                                                 "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                                                                 "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                                                                 "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))
subgenome_colors<-c("M1"="#7896BF", "M2"="#BF4342")
ggplot(subgenome_NA_Counts,aes(x= Subgenome, y= NA_Count))+
  geom_jitter(aes(color=Genome), height=0)+
  scale_color_manual(values = genome_colors)+
  theme_bw()+ ylab("Unaligned Exons Count")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/exonMissingness.bySubgenome.png", device="png",dpi=300, height = 6, width = 6)

filter(subgenome_NA_Counts, NA_Count >= 10000)
#M1: Zn, Zn4to1
#M2: TdKS, Zn4to1, ZdGigi, ZdMomo, Zn

#check that data is normally distributed to do two paired t test
shapiro.test(subgenome_NA_Counts$NA_Count)
#So not normally distributed...

pivot_wider(select(subgenome_NA_Counts, -name), values_from = NA_Count, names_from = Subgenome) %>% 
  mutate(difference = M1-M2) %>% 
  ggplot(aes(x=difference))+
  geom_histogram(bins = 35)+
  theme_bw()

#Here negative would mean that there's more NA's counted in M2 than M1?
pivot_wider(select(subgenome_NA_Counts, -name), values_from = NA_Count, names_from = Subgenome) %>% 
  mutate(difference = M1-M2) %>% select(difference) %>% t.test()
#Statistically significantly different 
#t = -13.492, df = 38, p-value = 4.534e-16
#alternative hypothesis: true mean is not equal to 0
#95 percent confidence interval:
#  -4706.678 -3478.553
#sample estimates:
#  mean of x 
#-4092.615  

mutate(subgenome_NA_Counts, NA_percent = (NA_Count/69269)*100) %>% 
  ggplot(aes(x= Subgenome, y= NA_percent))+
  geom_jitter(aes(color=Genome), height=0)+
  scale_color_manual(values = genome_colors)+
  theme_bw()+ ylab("Unaligned Exons Percentage")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/exonMissingnessPercent.bySubgenome.png", device="png",dpi=300, height = 6, width = 6)

mutate(subgenome_NA_Counts, NA_percent = (NA_Count/69269)*100) %>% 
  ggplot(aes(x= Subgenome, y= NA_percent))+
  geom_violin(aes(fill=Subgenome))+
  scale_fill_manual(values = subgenome_colors)+
  theme_bw()+ ylab("Unaligned Exons Percentage")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/exonMissingnessPercent.bySubgenome.violinplot.png", device="png",dpi=300, height = 6, width = 6)


#What about missingness by gene?
#How many genes have at least 1 uncalled exon?


#How often is there a gene that has some exons unaligned and some exons aligned?
genes_with_mixed_alignment<-tibble(Genome=NA, Subgenome=NA, Gene_ID=NA)
for(i in gene_id){
  df<-filter(full.fractionation.status, Gene_ID == i)
  for(t in tripsacinae_genome_IDs){
    for(m in c("M1","M2")){
      if(nrow(filter(df,any(is.na(get(paste0(t,".",m)))) & any(is.numeric(get(paste0(t,".",m)))))) > 0){
        genes_with_mixed_alignment<-add_row(genes_with_mixed_alignment, Genome = t, Subgenome = m, Gene_ID = i)
      }
    }
  }
}
genes_with_mixed_alignment<-na.omit(genes_with_mixed_alignment)
unique(genes_with_mixed_alignment$Gene_ID) %>% length() #9250

##Maize meeting figures##

#Create a figure that shows missing vs. alignment of exons by genome and subgenome
full.fractionation.status.summary<-full.fractionation.status %>% 
  select(c(contains("ID"),contains(".M"))) %>% 
  pivot_longer(cols = contains(".M"), names_to = "GenomeSubgenome", values_to = "Status") %>%
  mutate(Genome = str_split(GenomeSubgenome, "[.]", simplify = T)[,1] %>% factor(levels = c(c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                                                                            "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                                                                            "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                                                                            "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))),
         Subgenome = str_split(GenomeSubgenome, "[.]",simplify = T)[,2],
         Status_character = case_when(Status == 0 ~ "Retained", Status == 1 ~ "Deleted", is.na(Status) ~ "Unaligned"),
         Status_alignment = case_when(Status == 1 | Status == 0 ~ "Aligned", is.na(Status) ~ "Unaligned"))

full.fractionation.status.summary %>% group_by(Genome, Subgenome, Status_alignment) %>% count() %>%
  mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(x=Status_alignment, y=Percent, fill=Subgenome))+
  geom_violin()+
  #geom_jitter(height=0, aes(color=Genome))+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  scale_fill_manual(values=subgenome_colors)+
  scale_color_manual(values=genome_colors)+
  ylab("Percent of Exons")+xlab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/Exon_alignmentStatus.bySubgenome.png",device="png",dpi=300,width=6, height=8)

#Create a figure that shows missing vs. retained vs. deleted by genome and subgenome
full.fractionation.status.summary %>% group_by(Genome, Subgenome, Status_character) %>% count() %>%
  mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(x=Status_character, y=Percent, fill=Subgenome))+
  geom_violin()+
  #geom_jitter(height=0, aes(color=Genome))+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  scale_fill_manual(values=subgenome_colors)+
  scale_color_manual(values=genome_colors)+
  ylab("Percent of Exons")+xlab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/Exon_AllStatuses.bySubgenome.png",device="png",dpi=300,width=6, height=8)

##Let's try graphing the fractionation status as something other than barcharts

group_by(long.full.fractionation.status, Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  ggplot(aes(x=Status, y=Percent))+
  geom_violin()+
  geom_jitter(height = 0, width = 0.25, aes(color = Genome), size = 1, alpha = .75)+
  scale_color_manual(values = genome_colors)+
  theme_bw()+
  scale_x_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("Percentage")+xlab("")+
  ggtitle("Fractionation Status by Exon")+
  coord_flip()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.AllGenomes.violinandscatter.png",device="png",width = 7.5,height = 5)

group_by(long.full.fractionation.status, Status, Genome) %>% count() %>% mutate(Percent = (n/69269)*100) %>%
  filter(str_detect(Genome,"Zn",negate = T) & !Genome %in% c("ZdGigi","ZdMomo") ) %>%
  ggplot(aes(x=Status, y=Percent))+
  geom_violin()+
  geom_jitter(height = 0, width = 0.25, aes(color = Genome), size = 1, alpha = .75)+
  scale_color_manual(values = genome_colors)+
  theme_bw()+
  scale_x_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))+
  ylab("Percentage")+xlab("")+
  ggtitle("Fractionation Status by Exon")+
  coord_flip()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.NoZn.violinandscatter.png",device="png",width = 7.5,height = 5)

#Make table:
group_by(long.full.fractionation.status, Status, Genome) %>% 
  count() %>% 
  mutate(Percent = (n/69269)*100) %>% 
  pivot_wider(names_from = "Status", values_from = c("n", "Percent")) %>%
  write_tsv("Fractionation.Status.Table.tsv")

test_genome_colors<-c("ZmB73" = "#92ddb0",
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
####CONVERGENCE OF LOSS ACROSS GENES####

ref.gene<-full.fractionation.status$Gene_ID %>% unique()

df<-filter(full.fractionation.status, Gene_ID == ref.gene[1]) #filter fractionation status to get just 1 gene
exon.ct<-df$CDS_ID %>% unique() %>% length() #get the number of total exons

#for a gene, all exons are numeric, none are NA
df %>% select(ends_with(".M1")) %>% #select columns for fractionation status of M1
  lapply(function(x) any(is.na(x))) %>% #check to see if any column has even 1 NA
  as_tibble() %>% #convert to tibble
  pivot_longer(cols = everything()) %>% #make it an easier to work with format
  select(value) %>% pull() %>% as.numeric()%>% #change logical to numeric (FALSE = 0, TRUE = 1)
  sum() == 0 #if none of the genomes have an NA exon in the gene, then it should be 0, otherwise it will be greater than 0

#For a genome, this grabs all the exons that are "fractionated"
holder<-df %>% filter(get(paste0("TdFL",".M1"))==1) %>% 
  select(CDS_ID) %>% pull() %>% str_split(string=.,pattern="[.]",simplify = T) 

convergence<-tibble(Gene_ID=NA, ExonCt=NA,Genome=NA,M=NA, Loss_Pattern=NA)
gene_atleastOneNA<-tibble(Gene_ID=NA, ExonCt=NA,M=NA)
for(k in 1:length(ref.gene)){
  df<-filter(full.fractionation.status, Gene_ID == ref.gene[k])
  exon.ct<-df$CDS_ID %>% unique() %>% length()
  for(m in c(".M1",".M2")){
    if(df %>% select(ends_with(m)) %>% lapply(function(x) any(is.na(x))) %>% 
       as_tibble() %>% pivot_longer(cols=everything()) %>%
       select(value) %>% pull() %>% as.numeric() %>% sum() == 0){
      for(g in tripsacinae_genome_IDs){
        holder<-df %>% filter(get(paste0(g,m))==1)%>% select(CDS_ID) %>% pull() %>% 
          str_split(string=.,pattern="[.]",simplify = T)
        if(nrow(holder) > 0){
          convergence<-add_row(convergence, Gene_ID=ref.gene[k],
                               ExonCt =exon.ct, Genome=g,M=str_remove(m,"[.]"),
                               Loss_Pattern=paste(holder[,7],collapse = ":"))}
      }
    }else{
      gene_atleastOneNA<-add_row(gene_atleastOneNA, Gene_ID=ref.gene[k], ExonCt=exon.ct,M=str_remove(m,"[.]"))
    }
  }
}
#this will handle the NAs by not considering any gene that is NA in any genome

#Do it again but without Zn and Zd genomes
convergence.NoZnZd<-tibble(Gene_ID=NA, ExonCt=NA,Genome=NA,M=NA, Loss_Pattern=NA)
gene_atleastOneNA.NoZnZd<-tibble(Gene_ID=NA, ExonCt=NA,M=NA)
for(k in 1:length(ref.gene)){
  df<-filter(full.fractionation.status, Gene_ID == ref.gene[k])
  exon.ct<-df$CDS_ID %>% unique() %>% length()
  for(m in c(".M1",".M2")){
    if(df %>% select(ends_with(m)) %>% select(-c(starts_with("ZdGigi."),starts_with("ZdMomo."),contains("Zn"))) %>% lapply(function(x) any(is.na(x))) %>% 
       as_tibble() %>% pivot_longer(cols=everything()) %>%
       select(value) %>% pull() %>% as.numeric() %>% sum() == 0){
      for(g in c("TdKS","TdFL","ZdGigi_4to1","ZdMomo_4to1","ZhRIMHU001", "ZmB73", "ZmB97", "ZmCML103", "ZmCML228",
          "ZmCML247", "ZmCML277", "ZmCML322", "ZmCML333", "ZmCML52", "ZmCML69", "ZmHP301", "ZmIL14H", "ZmKi11",
          "ZmKi3", "ZmKy21", "ZmM162W", "ZmM37W", "ZmMo18W", "ZmMS71", "ZmNC350", "ZmNC358", "ZmOh43", "ZmOh7b",
          "ZmP39", "ZmTx303", "ZmTzi8",  "ZvTIL01", "ZvTIL11", "ZxTIL18", "ZxTIL25")){
        holder<-df %>% filter(get(paste0(g,m))==1)%>% select(CDS_ID) %>% pull() %>% 
          str_split(string=.,pattern="[.]",simplify = T)
        if(nrow(holder) > 0){
          convergence.NoZnZd<-add_row(convergence.NoZnZd, Gene_ID=ref.gene[k],
                               ExonCt =exon.ct, Genome=g,M=str_remove(m,"[.]"),
                               Loss_Pattern=paste(holder[,7],collapse = ":"))}
      }
    }else{
      gene_atleastOneNA.NoZnZd<-add_row(gene_atleastOneNA.NoZnZd, Gene_ID=ref.gene[k], ExonCt=exon.ct,M=str_remove(m,"[.]"))
    }
  }
}


#get rid of initializing NA rows
#convergence<-convergence[-1,]
convergence.NoZnZd<-convergence.NoZnZd[-1,]

#gene_atleastOneNA<-gene_atleastOneNA[-1,]
gene_atleastOneNA.NoZnZd<-gene_atleastOneNA.NoZnZd[-1,]

#How much of each gene was lost?
#filter(convergence, Gene_ID == convergence$Gene_ID[1]) %>% mutate(NumExonsInPattern = str_count(Loss_Pattern, pattern = ":")+1,
#                                                                  PercentGeneLost = (NumExonsInPattern/ExonCt)*100)

#convergence<-convergence %>% mutate(NumExonsInPattern = str_count(Loss_Pattern, pattern = ":")+1,
#                                    PercentGeneLost = (NumExonsInPattern/ExonCt)*100)
convergence.NoZnZd<-convergence.NoZnZd %>% mutate(NumExonsInPattern = str_count(Loss_Pattern, pattern = ":")+1,
                                    PercentGeneLost = (NumExonsInPattern/ExonCt)*100)

#convergence %>% group_by(Genome, M) %>% 
#  reframe(max.PerGeneLost=max(PercentGeneLost), min.PerGeneLost=min(PercentGeneLost),
#          median.PerGeneLost=median(PercentGeneLost),mean.PerGeneLost=mean(PercentGeneLost)) %>% print(n=70)

#convergence$M<-convergence$M %>% factor()
#convergence$Genome<-convergence$Genome %>% factor()
convergence.NoZnZd$M<-convergence.NoZnZd$M %>% factor()
convergence.NoZnZd$Genome<-convergence.NoZnZd$Genome %>% factor()
#PerecentGeneLost.model<-aov(PercentGeneLost ~ Genome + M, convergence)
#PerGeneLost.M.contrast<-multcomp::glht(PerecentGeneLost.model, linfct = multcomp::mcp(M="Tukey"))
#PerGeneLost.M.letters<-multcomp::cld(PerGeneLost.M.contrast)
#PerGeneLost.Genome.contrast<-multcomp::glht(PerecentGeneLost.model, linfct = multcomp::mcp(Genome="Tukey"))
#PerGeneLost.Genome.letters<-multcomp::cld(PerGeneLost.Genome.contrast)
#PerGeneLost.Genome.letters<-TukeyHSD(PerecentGeneLost.model, which = "Genome")
#PerGeneLost.Genome.letters$Genome %>% as_tibble() %>% add_column(contrast=PerGeneLost.Genome.letters$Genome %>% dimnames() %>% .[[1]])%>%
#  filter(`p adj` < 0.05) %>% print(n=40)
#Everything vs. TdFL passes the p adj fileter, nothing else

#ggplot(convergence, aes(fill=M, x= PercentGeneLost, y=Genome))+
#  geom_boxplot()+
#  theme_minimal()

#ggplot(convergence, aes(fill=M, x=PercentGeneLost))+
#  facet_wrap(facets = vars(Genome))+
#  geom_histogram()

if(!"Gene_ID" %in% colnames(ref_Sb313.cds)){#if Gene_ID is not a column in the ref_Sb313.cds object
  #add it as a column
  ref_Sb313.cds<-ref_Sb313.cds %>% mutate(Gene_ID = str_split(ID,";",simplify = T)[,2] %>% str_remove_all(pattern = "Parent="))
}
#convergence<-convergence %>% 
#  mutate(RefChr = case_when(Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "01")[,7]) ~ "Chr01",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "02")[,7]) ~ "Chr02",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "03")[,7]) ~ "Chr03",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "04")[,7]) ~ "Chr04",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "05")[,7]) ~ "Chr05",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "06")[,7]) ~ "Chr06",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "07")[,7]) ~ "Chr07",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "08")[,7]) ~ "Chr08",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "09")[,7]) ~ "Chr09",
#                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "10")[,7]) ~ "Chr10"))
convergence.NoZnZd<-convergence.NoZnZd %>% 
  mutate(RefChr = case_when(Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "01")[,7]) ~ "Chr01",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "02")[,7]) ~ "Chr02",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "03")[,7]) ~ "Chr03",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "04")[,7]) ~ "Chr04",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "05")[,7]) ~ "Chr05",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "06")[,7]) ~ "Chr06",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "07")[,7]) ~ "Chr07",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "08")[,7]) ~ "Chr08",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "09")[,7]) ~ "Chr09",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "10")[,7]) ~ "Chr10"))
#convergence$Genome<-convergence$Genome %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                             "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25",
#                                                             "ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))
convergence.NoZnZd$Genome<-convergence.NoZnZd$Genome %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                                             "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25",
                                                             "ZhRIMHU001","ZdGigi_4to1","ZdMomo_4to1" ,"TdFL","TdKS"))

library(ggridges)
#filter(convergence, !is.na(RefChr)) %>%
#  ggplot(aes(x=PercentGeneLost))+
#  geom_density_ridges(aes(y=Genome,fill=Genome))+
#  theme_minimal()+
#  scale_fill_manual(values = genome_colors)+
#  facet_wrap(facets = vars(RefChr))+
#  xlab("Percent Gene Lost")+ylab("")+
#  ggtitle("Percent Gene Lost by Sb Ref Chromosome")
#scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained"))
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/PerGeneLost.ByRefChr.AllGenomes.png", device = "png",dpi=500,height=7, width=8)

#filter(convergence, !is.na(RefChr)) %>%
#  ggplot(aes(x=PercentGeneLost))+
#  geom_density_ridges(aes(y=Genome,fill=Genome))+
#  theme_minimal()+
#  scale_fill_manual(values = genome_colors)+
#  xlab("Percent Gene Lost")+ylab("")+
#  ggtitle("Percent Gene Lost")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/PerGeneLost.AllGenomes.png", device = "png",dpi=500,height=7, width=8)

#filter(convergence, !is.na(RefChr)) %>%
#  ggplot(aes(x=PercentGeneLost))+
#  geom_density_ridges(aes(y=Genome,fill=Genome))+
#  theme_bw()+
#  facet_wrap("M")+
#  scale_fill_manual(values = genome_colors)+
#  xlab("Percent Gene Lost")+ylab("")+
#  ggtitle("Percent Gene Lost")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/PerGeneLost.AllGenomes.BySubgenome.png",device = "png",dpi=500,height=7, width=8)

filter(convergence.NoZnZd, !is.na(RefChr)) %>%
  ggplot(aes(x=PercentGeneLost))+
  geom_density_ridges(aes(y=Genome,fill=Genome))+
  theme_bw()+
  facet_wrap("M")+
  scale_fill_manual(values = genome_colors)+
  xlab("Percent Gene Lost")+ylab("")+
  ggtitle("Percent Gene Lost")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/PerGeneLost.NoZnZd.BySubgenome.png",device = "png",dpi=500,height=7, width=8)

filter(convergence.NoZnZd, !is.na(RefChr)) %>% group_by(M, Genome) %>% reframe(perLoss.mean = mean(PercentGeneLost, na.rm = T)) %>%
  ggplot(aes(x=M, y=perLoss.mean, fill=M))+
  geom_violin()+
  #geom_text(aes(label = Genome))+
  scale_fill_manual(values = subgenome_colors)+
  theme_bw()+ylab("Mean Percentage of Gene Lost")+xlab("Subgenome")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/PerGeneLost.NoZnZd.BySubgenome.mean.png",device="png",dpi=300,height = 5,width = 5)

#filter(convergence, !is.na(RefChr)) %>% group_by(M, Genome) %>% reframe(perLoss.mean = mean(PercentGeneLost, na.rm = T)) %>%
#  ggplot(aes(x=M, y=perLoss.mean, fill=M))+
#  geom_violin()+
#  #geom_text(aes(label = Genome))+
#  scale_fill_manual(values = subgenome_colors)+
#  theme_bw()+ylab("Mean Percentage of Gene Lost")+xlab("Subgenome")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/PerGeneLost.AllGenomes.BySubgenome.mean.png",device="png",dpi=300,height = 5,width = 5)

filter(convergence.NoZnZd, !is.na(RefChr)) %>% group_by(M, Genome) %>% reframe(perLoss.mean = mean(PercentGeneLost, na.rm = T), perLoss.median = median(PercentGeneLost, na.rm = T)) %>% view()
filter(convergence.NoZnZd, !is.na(RefChr)) %>% group_by(M, Genome) %>% reframe(perLoss.min = min(PercentGeneLost, na.rm=T))%>% view()
#aov(PercentGeneLost ~ M+Genome, convergence) %>% summary()
#Df    Sum Sq Mean Sq  F value  Pr(>F)    
#M                1  10351339 10351339 8084.121  < 2e-16 ***
#Genome          38    129310     3403    2.658 1.28e-07 ***
#  Residuals   397262 508675415     1280     
#TukeyHSD(aov(PercentGeneLost ~ M+Genome, convergence))

#TukeyHSD(aov(PercentGeneLost ~ M+Genome, convergence), which = "Genome") %>% .$Genome %>% view()
#Only significant differences are between Td genomes and Zea genomes
#TukeyHSD(aov(PercentGeneLost ~ RefChr+M+Genome, convergence), which = "RefChr") %>% .$RefChr %>% view()

aov(PercentGeneLost ~ M+Genome, convergence.NoZnZd) %>% summary()
#Df    Sum Sq  Mean Sq   F value Pr(>F)    
#M                1  14576761 14576761 11528.30 <2e-16 ***
#Genome          34    200761     5905     4.67 <2e-16 ***
#  Residuals   516235 652744586     1264 

#convergence %>% group_by(RefChr) %>% summarize(avg=mean(PercentGeneLost,na.rm=T), median=median(PercentGeneLost,na.rm=T)) 

#how many loss patterns per gene?
#filter(convergence, Gene_ID == convergence$Gene_ID[1]) %>% group_by(M) %>% reframe(uniqLossPattern = unique(Loss_Pattern))

#convergence.summary<-tibble(Gene_ID=NA, M=NA, Number_Unique_Loss_Patterns=NA)
#for(i in unique(convergence$Gene_ID)){
#  df<-filter(convergence, Gene_ID == i) %>% group_by(M) %>% reframe(uniqLossPattern = unique(Loss_Pattern))
#  M1.uniq.pattern<-filter(df, M == "M1") %>% nrow()
#  M2.uniq.pattern<-filter(df, M == "M2") %>% nrow()
#  if(M1.uniq.pattern > 0 & !is.na(M1.uniq.pattern)){
#    convergence.summary<-add_row(convergence.summary, Gene_ID=i, M="M1",Number_Unique_Loss_Patterns=M1.uniq.pattern)
#  }
#  if(M2.uniq.pattern > 0 & !is.na(M2.uniq.pattern)){
#    convergence.summary<-add_row(convergence.summary, Gene_ID=i, M="M2",Number_Unique_Loss_Patterns=M2.uniq.pattern)
#  }
#}

convergence.summary.NoZdZn<-tibble(Gene_ID=NA, M=NA, Number_Unique_Loss_Patterns=NA)
for(i in unique(convergence.NoZnZd$Gene_ID)){
  df<-filter(convergence.NoZnZd, Gene_ID == i) %>% group_by(M) %>% reframe(uniqLossPattern = unique(Loss_Pattern))
  M1.uniq.pattern<-filter(df, M == "M1") %>% nrow()
  M2.uniq.pattern<-filter(df, M == "M2") %>% nrow()
  if(M1.uniq.pattern > 0 & !is.na(M1.uniq.pattern)){
    convergence.summary.NoZdZn<-add_row(convergence.summary.NoZdZn, Gene_ID=i, M="M1",Number_Unique_Loss_Patterns=M1.uniq.pattern)
  }
  if(M2.uniq.pattern > 0 & !is.na(M2.uniq.pattern)){
    convergence.summary.NoZdZn<-add_row(convergence.summary.NoZdZn, Gene_ID=i, M="M2",Number_Unique_Loss_Patterns=M2.uniq.pattern)
  }
}

#convergence.summary<-convergence.summary[-1,]
convergence.summary.NoZdZn<-convergence.summary.NoZdZn[-1,]

#ggplot(convergence.summary, aes(x=Number_Unique_Loss_Patterns, fill=M))+
#  geom_histogram(binwidth = 1)+
#  theme_minimal()+
#  scale_fill_manual(values = subgenome_colors)+
#  xlab("Number of ways a gene was lost")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/SummaryUniqPatternsLoss.ByGene.AllGenomes.histogram.png", device="png",dpi=300,width = 4, height = 4.5)

#ggplot(convergence.summary, aes(y=Number_Unique_Loss_Patterns, x=M, fill=M))+
#  geom_violin()+
#  ylab("")+xlab("")+ggtitle("Number of ways a gene was lost")+
#  theme_minimal()+
#  scale_fill_manual(values = subgenome_colors)+
#  theme(legend.position = "none")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/SummaryUniqPatternsLoss.ByGene.AllGenomes.violin.png", device="png",dpi=300,width = 4, height = 4.5)

#convergence.summary %>% group_by(M) %>% reframe(mean.LossPatterns = mean(Number_Unique_Loss_Patterns),
#                                                median.LossPatterns = median(Number_Unique_Loss_Patterns),
#                                                sd.LossPatterns = sd(Number_Unique_Loss_Patterns),
#                                                min.LossPatterns = min(Number_Unique_Loss_Patterns),
#                                                max.LossPatterns = max(Number_Unique_Loss_Patterns))


convergence.summary.NoZdZn %>% group_by(M) %>% reframe(mean.LossPatterns = mean(Number_Unique_Loss_Patterns),
                                                       median.LossPatterns = median(Number_Unique_Loss_Patterns),
                                                       sd.LossPatterns = sd(Number_Unique_Loss_Patterns),
                                                       min.LossPatterns = min(Number_Unique_Loss_Patterns),
                                                       max.LossPatterns = max(Number_Unique_Loss_Patterns))


#Is there a significant difference in the number of unique loss patterns identified for M1 vs. M2
#This will only be comparing unique loss patterns if both M1 and M2 have a loss pattern (if retained or unaligned in 1 subgenome, it won't be included here)
#convergence.summary %>% pivot_wider(id_cols = "Gene_ID", names_from = "M", values_from = "Number_Unique_Loss_Patterns")%>%
#  na.omit()%>% #loses 6K (8649 --> 2423)
  #pivot_longer(cols=starts_with("M"),names_to = "M",values_to = "Number_Unique_Loss_Patterns")%>%
  #not testing if these are normally distributed
#  mutate(difference = M1 - M2)%>%
#  select(difference) %>%
#  t.test()
#t = 0.94677, df = 2422, p-value = 0.3438
#alternative hypothesis: true mean is not equal to 0
#95 percent confidence interval:
#  -0.02652556  0.07605094
#sample estimates:
#  mean of x 
#0.02476269


#Originally suggested that there are more unique loss patterns in M1 than M2 (significantly), but that goes away when including other genomes and not really (difference is 0.16)
#See /work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/SummaryUniqPatternsLoss.ByGene.AllGenomes.violin.png

convergence.summary.NoZdZn %>% pivot_wider(id_cols = "Gene_ID", names_from = "M", values_from = "Number_Unique_Loss_Patterns")%>%
  na.omit()%>% # 11119 --> 7299
  #pivot_longer(cols=starts_with("M"),names_to = "M",values_to = "Number_Unique_Loss_Patterns")%>%
  #not testing if these are normally distributed
  mutate(difference = M1 - M2)%>%
  select(difference) %>%
  t.test()
#t = 1.3542, df = 5514, p-value = 0.1757
#alternative hypothesis: true mean is not equal to 0
#95 percent confidence interval:
#  -0.01022624  0.05591980
#sample estimates:
#  mean of x 
#0.02284678 

#Same pattern as before, but the magnitude of difference is smaller when excluding Zn and Zd

ggplot(convergence.summary.NoZdZn, aes(y=Number_Unique_Loss_Patterns, x=M, fill=M))+
  geom_violin()+
  ylab("")+xlab("")+ggtitle("Number of ways a gene was lost")+
  theme_minimal()+
  scale_fill_manual(values = subgenome_colors)+
  theme(legend.position = "none")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/SummaryUniqPatternsLoss.ByGene.NoZnZd.violin.png", device="png",dpi=300,width = 4, height = 4.5)


#NumLossPattterns.model<-aov(Number_Unique_Loss_Patterns ~ M, convergence.summary)
#summary(NumLossPattterns.model)
#               Df Sum Sq Mean Sq F value Pr(>F)    
#M               1      0  0.2996   0.282  0.596
#Residuals   11233  11955  1.0643 
#TukeyHSD(NumLossPattterns.model)
#diff        lwr        upr p adj
#M2-M1 -0.01062889 -0.04989415 0.02863638 0.5957006
NumLossPattterns.NoZdZn.model<-aov(Number_Unique_Loss_Patterns ~ M, convergence.summary.NoZdZn)
summary(NumLossPattterns.NoZdZn.model)
#Df Sum Sq Mean Sq F value Pr(>F)
#M               1      1  1.4093   1.416  0.234
#Residuals   16161  16090  0.9956
TukeyHSD(NumLossPattterns.NoZdZn.model)
#Same finding as the t-tests above really
#diff         lwr       upr     p adj
#M2-M1 -0.01881586 -0.04981484 0.01218311 0.2341603

#Pick out a couple of genes that could be good to highlight
convergence.summary.NoZdZn %>% filter(Number_Unique_Loss_Patterns > 10)
filter(convergence.NoZnZd, Gene_ID == "Sobic.001G011700.1.v3.1" & M == "M2") %>% 
  write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/exampleGeneForConvergence.tsv")

#What about genes with at least one genome having at least 1 unaligned exon
#gene_atleastOneNA %>% group_by(M) %>% count() #M1=4304 and M2=7507
gene_atleastOneNA.NoZnZd %>% group_by(M) %>% count() #M1=1818 and M2=4439

#how often do different genomes lose the same gene in different ways?
#note to self, none of the loss_patterns are NA, none are ""
#let's make the loss_pattern sharing a proportion
#if they share the exact same pattern = 1
#if they share some exons, number exons shared / total exon cnt
#if no sharing, 0
#convergence.subset
#write_tsv(convergence[1:3500,], file = "/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_test/convergence.tsv")
write_tsv(convergence.NoZnZd, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_test/full.convergence.NoZnZd.tsv")
write_tsv(gene_atleastOneNA.NoZnZd, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_test/gene_atleastOneNA.NoZnZd.tsv")

#convergence.shared<-tibble(Gene_ID=NA, genome1=NA, genome2=NA, M=NA, Loss_Pattern_Sharing_Per=NA,Loss_Pattern_Sharing_Abs=NA)

#test<-filter(convergence, Gene_ID == "Sobic.001G004700.1.v3.1")

# should the exact match percent sharing also be divided by the total number of exons in the gene?

#for(g in 1:38){ #for each genome
#  for(counter in (g+1):39){ #paired with each other genome
#    for(b in c("M1","M2")){ #and for each subgenome
#      if(nrow(filter(gene_atleastOneNA, Gene_ID == test$Gene_ID[1] & M == b))==0){ #if the gene x subgenome combo isn't in the list of genes with an NA
#        pattern1<-filter(test, Genome == tripsacinae_genome_IDs[g] & M == b) %>% select(Loss_Pattern) %>% pull() 
#        pattern2<-filter(test, Genome == tripsacinae_genome_IDs[counter] & M == b) %>% select(Loss_Pattern) %>% pull() 
#        if(is_empty(pattern1) | is_empty(pattern2)){ #if one of the patterns is empty (AKA the gene is completely retained in one of the genomes)
#          convergence.shared<-add_row(convergence.shared, Gene_ID = test$Gene_ID[1], genome1=tripsacinae_genome_IDs[g], genome2=tripsacinae_genome_IDs[counter], 
#                                                         M=b, Loss_Pattern_Sharing_Per=NA,Loss_Pattern_Sharing_Abs=NA) #NAs indicate the gene is retained in one of the genomes
#        }else{
#          if(pattern1 == pattern2){ #if the loss patterns are exactly the same
#            convergence.shared<-add_row(convergence.shared, Gene_ID = test$Gene_ID[1], genome1=tripsacinae_genome_IDs[g], genome2=tripsacinae_genome_IDs[counter], 
#                                                              M=b, Loss_Pattern_Sharing_Per=(str_count(pattern1, pattern = ":")+1)/pull(select(filter(test, Genome == tripsacinae_genome_IDs[g] & M == b),ExonCt)),Loss_Pattern_Sharing_Abs=1)
#            #Loss sharing percent 1 == complete sharing, Loss pattern sharing ABS == binary for complete sharing (1 exact share, 0 not exact share)
#          }else{ #if there's a pattern for both genomes, but they're not exactly the same
#            if(length(str_split(pattern1, ":",simplify = T)) >= length(str_split(pattern2,":",simplify = T))){
#              longest.pattern<-pattern1
#              shortest.pattern<-pattern2
#              #figure out which one is the shortest and longest pattern
#            }else{
#              longest.pattern<-pattern2
#              shortest.pattern<-pattern1
#            }
#            shortest.pattern<-str_split(shortest.pattern,":", simplify=T)
#            convergence.shared<-add_row(convergence.shared, Gene_ID = test$Gene_ID[1], genome1=tripsacinae_genome_IDs[g], genome2=tripsacinae_genome_IDs[counter],
#                                        M=b, Loss_Pattern_Sharing_Per=sum(str_count(longest.pattern, str_c(shortest.pattern, "\\b")))/pull(select(filter(test, Genome == tripsacinae_genome_IDs[g] & M == b),ExonCt)),
#                                        Loss_Pattern_Sharing_Abs=0) 
#            #Loss_pattern_Sharing_Per is counting the number of entries in the shortest pattern that can be found in the longest pattern divided by the total exon count
#            #The "\\b" part keeps it from double counting (ex: 1 being counted for 1, 10, 11, etc.)
#          }
#        }
#      }
#    }
#    counter<-counter+1
#  }
#}
#convergence.shared<-convergence.shared[-1,]
#Code from Arun to parallalize convergence.shared loop is in convergence_loop.R
####Arun ran it separately and created the full object here: convergence_shared.v3.tsv ####
#convergence.shared<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_test/convergence_shared.v3.tsv")

#select(convergence.shared, contains("genome")) %>% unique()

#convergence.shared$genome1<-convergence.shared$genome1 %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                                             "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,#"ZxTIL25",
#                                                                             "ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))
#convergence.shared$genome2<-convergence.shared$genome2 %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                                             "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL"))

#convergence.shared.stats<-convergence.shared %>% group_by(M, genome1, genome2) %>% 
#  reframe(mean.LossSharePer = mean(Loss_Pattern_Sharing_Per,na.rm = T), #how many of the same exons are lost in genome pair/total exons?
#          mean.LossShareAbs = mean(Loss_Pattern_Sharing_Abs, na.rm=T), #are they sharing the exact same exons lost
#          median.LossSharePer = median(Loss_Pattern_Sharing_Per,na.rm = T),
#          median.LossShareAbs = median(Loss_Pattern_Sharing_Abs, na.rm=T),
#          cnt.NA = sum(is.na(Loss_Pattern_Sharing_Abs)),
#          sum.LossShareAbs = sum(Loss_Pattern_Sharing_Abs, na.rm=T))

#ggplot(data = convergence.shared.stats, aes(x=genome1, y=genome2, fill=mean.LossShareAbs))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
#  #geom_text(aes(label=sprintf("%.2f",mean.LossShareAbs)), size=2, color="darkgrey")+
#  ggtitle("Average Absolute Sharing of Gene Loss Pattern")+#xlab("")+ylab("")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+
#  scale_fill_viridis_c(option = "magma", direction = 1, name = "mean")

#convergence.shared.stats <- convergence.shared.stats %>% mutate(genome1_graph = case_when(genome1 %in% c("ZhRIMHU001") & 
#                                                                                            genome2 %in% c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                                                                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                                                           "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                                                           "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18") ~ genome2,
#                                                                                          genome1 %in% c("ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69") & 
#                                                                                            genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                                                                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome2,
#                                                                                          genome1 == "ZmHP301" & genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome2,
#                                                                                          genome1 %in% c("ZmM37W","ZmMo18W") & genome2 %in% c("ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome2, 
#                                                                                          genome1 %in% c("ZmKi11","ZmKi3") & 
#                                                                                            genome2 %in% c("ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39","ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome2,
#                                                                                          genome1 %in% c("ZmNC350","ZmNC358") & genome2 %in% c("ZmOh43","ZmOh7b","ZmP39","ZmTx303") ~ genome2,
#                                                                                          genome1 == "ZvTIL01" & genome2 == "ZvTIL11" ~ genome2,
#                                                                                          .default = genome1
#),
#genome2_graph = case_when(genome1 %in% c("ZhRIMHU001", "ZdMomo","ZdGigi") & 
#                            genome2 %in% c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                           "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                           "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18") ~ genome1,
#                          #genome2 == "ZhRIMHU001" & genome1 %in% c("ZdMomo","ZdGigi") ~ genome1,
#                          genome1 == "ZnPI615697" & genome2 %in% c("ZvTIL01","ZvTIL11","ZxTIL18") ~ genome1, 
#                          genome1 %in% c("ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69") & 
#                            genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome1,
#                          genome1 == "ZmHP301" & genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome1,
#                          genome1 %in% c("ZmM37W","ZmMo18W") & genome2 %in% c("ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome1, 
#                          genome1 %in% c("ZmKi11","ZmKi3") & 
#                            genome2 %in% c("ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39","ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome1,
#                          genome1 %in% c("ZmNC350","ZmNC358") & genome2 %in% c("ZmOh43","ZmOh7b","ZmP39","ZmTx303") ~ genome1,
#                          genome1 == "ZvTIL01" & genome2 == "ZvTIL11" ~ genome1,
#                          .default = genome2
#)
#)
#convergence.shared.stats$genome1_graph<-convergence.shared.stats$genome1_graph %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                                                                     "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                                                     "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                                                     "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,#"ZxTIL25",
#                                                                                                     "ZhRIMHU001" ,"ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))
#convergence.shared.stats$genome2_graph<-convergence.shared.stats$genome2_graph %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                                                                     "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                                                     "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                                                     "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))
#
#convergence.shared.stats %>% filter(genome1 == "TdFL") %>% arrange(-sum.LossShareAbs)
#convergence.shared.stats %>% filter(genome1 == "TdFL") %>% ungroup() %>% reframe(mean=mean(sum.LossShareAbs))
#convergence.shared.stats %>% filter(genome1 != "TdFL") %>% arrange(sum.LossShareAbs)
#convergence.shared.stats %>% filter(genome1 != "TdFL") %>% ungroup() %>% reframe(mean=mean(sum.LossShareAbs))
#convergence.shared.stats %>% filter(!str_detect(genome1,"Zm") & !str_detect(genome2,"Zm") & genome1 != "TdFL") %>% ungroup() %>% reframe(mean=mean(sum.LossShareAbs))

#group_by(convergence.shared.stats, M) %>% reframe(cnt.NA.mean=mean(cnt.NA, na.rm=T))
#difference in average differential retnetion between M1 and M2
# M1 average = 1925, M2 average = 1193
##plots with all genomes
#ggplot(data = convergence.shared.stats, aes(x=genome1_graph, y=genome2_graph, fill=mean.LossShareAbs))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
#  #geom_text(aes(label=sprintf("%.2f",mean.LossShareAbs)), size=2, color="darkgrey")+
#  ggtitle("Average Absolute Sharing of Gene Loss Pattern")+xlab("")+ylab("")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+
#  scale_fill_viridis_c(option = "magma", direction = 1, name = "mean")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.allgenomes.avgAbsLoss.png", device="png",dpi=300,width = 15, height=8)
#ggplot(data = convergence.shared.stats, aes(x=genome1_graph, y=genome2_graph, fill=sum.LossShareAbs))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
##  geom_text(aes(label=sum.LossShareAbs), size=2, color="darkgrey")+
#  ggtitle("Total Absolute Sharing of Gene Loss Pattern")+xlab("")+ylab("")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+
#  scale_fill_viridis_c(option = "magma", direction = 1, name = "count")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.allgenomes.cntAbsLoss.png", device="png",dpi=300,width = 15, height=8)

#ggplot(data = convergence.shared.stats, aes(x=genome1_graph, y=genome2_graph, fill=mean.LossSharePer))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
##  geom_text(aes(label=sprintf("%.2f",mean.LossSharePer)), size=2, color="darkgrey")+
#  ggtitle("Average Percent Sharing of Gene Loss Pattern")+xlab("")+ylab("")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+
#  scale_fill_viridis_c(option = "magma", direction = 1, name = "mean")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.allgenomes.avgPerLoss.png", device="png",dpi=300,width = 15, height=8)

#ggplot(data=convergence.shared.stats, aes(x=genome1_graph, y=genome2_graph, fill=cnt.NA))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
##  geom_text(aes(label=cnt.NA), size=2, color="darkgrey")+
#  ggtitle("Total in retained vs. some loss, pairwise")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+
#  scale_fill_viridis_c(option = "magma", direction = 1, name = "count")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.allgenomes.cntNA.png", device="png",dpi=300,width = 15, height=8)

##Maize meeting plots, turn these to percentages##
#12157 is the number of reference genes
#mutate(convergence.shared.stats, total.NA.percent = (cnt.NA/12157)*100, total.LossShareAbs.percent = (sum.LossShareAbs/12157)*100)%>%
#  #filter(!genome1 %in% c("ZnPI615697","ZdGigi","ZdMomo","TdFL") &!genome2 %in% c("ZnPI615697","ZdGigi","ZdMomo")) %>%
#  ggplot(aes(x=genome1_graph, y=genome2_graph, fill=total.NA.percent))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
#  scale_fill_viridis_c(option="magma",direction=1, name = "Percent")+
#  ggtitle("Total percent of genes showing completely retained in one genome vs. some loss")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+xlab("")+ylab("")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.PerTotalRetainedVSSomeLoss.png", device = "png",dpi=300, width = 15, height=8)

#convergence.shared.stats %>% filter(genome1 == "TdFL") %>% reframe(avg=mean(cnt.NA), min=min(cnt.NA),max=max(cnt.NA))
##the range of completely retained vs. some loss counts for TdFL vs. anything
##avg = 1938 genes, min = 1524, max = 2342 genes
#convergence.shared.stats %>% filter(genome1 == "TdFL") %>% group_by(M)%>% reframe(avg=mean(cnt.NA), min=min(cnt.NA),max=max(cnt.NA))
#M1: average=2322 genes, min=2299, max=2342
#M2: average=1554 genes, min=1524, max=1579
#convergence.shared.stats %>% group_by(M) %>% reframe(avg=mean(cnt.NA), min=min(cnt.NA),max=max(cnt.NA))

#convergence.shared <-convergence.shared %>% mutate(genome.paired = str_c(genome1,genome2, sep=":"))

#LossAbs.share.model<-aov(Loss_Pattern_Sharing_Abs ~ genome.paired, convergence.shared)
#summary(LossAbs.share.model)
#Df Sum Sq Mean Sq F value Pr(>F)    
#
#LossAbs.shareZhmvx.model<-aov(Loss_Pattern_Sharing_Abs ~ genome.paired, filter(convergence.shared, !genome1 %in% c("ZnPI615697","ZdGigi","ZdMomo","TdFL") &!genome2 %in% c("ZnPI615697","ZdGigi","ZdMomo")))
#summary(LossAbs.shareZhmvx.model)
#Df Sum Sq Mean Sq F value Pr(>F)    
#

#LossPer.share.model<-aov(Loss_Pattern_Sharing_Per ~ genome.paired, convergence.shared)
#summary(LossPer.share.model)
#Df Sum Sq Mean Sq F value Pr(>F)    
#
#LossPer.shareZhmvx.model.<-aov(Loss_Pattern_Sharing_Per ~ genome.paired, filter(convergence.shared, !genome1 %in% c("ZnPI615697","ZdGigi","ZdMomo","TdFL") &!genome2 %in% c("ZnPI615697","ZdGigi","ZdMomo")))
#summary(LossPer.shareZhmvx.model.)
#Df Sum Sq Mean Sq F value Pr(>F)    
#

#add in the information about each gene in each genome too
#convergence.shared<-
#  left_join(x=convergence.shared, 
#          y=convergence.NoZnZd,
#          by=c("Gene_ID","M", "genome1" = "Genome"))
#convergence.shared<-
#  left_join(x=convergence.shared, 
#            y=convergence.NoZnZd,
#            by=c("Gene_ID","M", "genome2" = "Genome"),
#            suffix = c(".genome1",".genome2"))

#colnames(convergence.shared)

#remove redundant columns and rename
#convergence.shared<- convergence.shared %>% select(-c("ExonCt.genome2","RefChr.genome2"))
#colnames(convergence.shared)[7]<-"ExonCt"
#colnames(convergence.shared)[11]<-"RefChr"

#Trying to figure out what exactly I want to look at
#filter(convergence.shared, Loss_Pattern_Sharing_Abs != 1 & genome1 == "TdFL") %>% view()

#What about the difference in percent loss? 

#convergence.shared<-mutate(convergence.shared, PercentGeneLost.diff = case_when(
#                          is.na(PercentGeneLost.genome1) != T & is.na(PercentGeneLost.genome2) != T ~ 
#                            PercentGeneLost.genome1 - PercentGeneLost.genome2,
#                          is.na(PercentGeneLost.genome1) == T & is.na(PercentGeneLost.genome2) != T ~ 
#                            0 - PercentGeneLost.genome2,
#                          is.na(PercentGeneLost.genome1) != T & is.na(PercentGeneLost.genome2) == T ~ 
#                            PercentGeneLost.genome1 - 0,
#                          is.na(PercentGeneLost.genome1) == T & is.na(PercentGeneLost.genome2) == T ~ 
#                            0)) 
#convergence.shared %>%
#  group_by(M, genome1,genome2) %>%
#  filter(PercentGeneLost.diff != 0) %>%
#  reframe(avg.PercentGeneLost.diff = mean(PercentGeneLost.diff, na.rm =T)) %>%
#  ggplot(aes(x=genome1, y=genome2, fill = avg.PercentGeneLost.diff))+
#  geom_tile()+
#  scale_fill_viridis_c(option="magma",direction=1, name = "Average")

#convergence.shared.stats<-convergence.shared %>% group_by(M, genome1, genome2) %>% 
#  reframe(mean.LossSharePer = mean(Loss_Pattern_Sharing_Per,na.rm = T), #how many of the same exons are lost in genome pair/total exons?
#          mean.LossShareAbs = mean(Loss_Pattern_Sharing_Abs, na.rm=T), #are they sharing the exact same exons lost
#          median.LossSharePer = median(Loss_Pattern_Sharing_Per,na.rm = T),
#          median.LossShareAbs = median(Loss_Pattern_Sharing_Abs, na.rm=T),
#          cnt.NA = sum(is.na(Loss_Pattern_Sharing_Abs)),
#          sum.LossShareAbs = sum(Loss_Pattern_Sharing_Abs, na.rm=T),
#          mean.allgenes.PercentGeneLost.diff = mean(PercentGeneLost.diff, na.rm=T))


#convergence.shared.stats <- convergence.shared.stats %>% mutate(genome1_graph = case_when(genome1 %in% c("ZhRIMHU001") & 
#                                                                                            genome2 %in% c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                                                                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                                                           "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                                                           "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18") ~ genome2,
#                                                                                          genome1 %in% c("ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69") & 
#                                                                                            genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                                                                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome2,
#                                                                                          genome1 == "ZmHP301" & genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome2,
#                                                                                          genome1 %in% c("ZmM37W","ZmMo18W") & genome2 %in% c("ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome2, 
#                                                                                          genome1 %in% c("ZmKi11","ZmKi3") & 
#                                                                                            genome2 %in% c("ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39","ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome2,
#                                                                                          genome1 %in% c("ZmNC350","ZmNC358") & genome2 %in% c("ZmOh43","ZmOh7b","ZmP39","ZmTx303") ~ genome2,
#                                                                                          genome1 == "ZvTIL01" & genome2 == "ZvTIL11" ~ genome2,
#                                                                                          .default = genome1
#),
#genome2_graph = case_when(genome1 %in% c("ZhRIMHU001", "ZdMomo","ZdGigi") & 
#                            genome2 %in% c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                           "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                           "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18") ~ genome1,
#                          #genome2 == "ZhRIMHU001" & genome1 %in% c("ZdMomo","ZdGigi") ~ genome1,
#                          genome1 == "ZnPI615697" & genome2 %in% c("ZvTIL01","ZvTIL11","ZxTIL18") ~ genome1, 
#                          genome1 %in% c("ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69") & 
#                            genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39",
#                                           "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome1,
#                          genome1 == "ZmHP301" & genome2 %in% c("ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome1,
#                          genome1 %in% c("ZmM37W","ZmMo18W") & genome2 %in% c("ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39") ~ genome1, 
#                          genome1 %in% c("ZmKi11","ZmKi3") & 
#                            genome2 %in% c("ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39","ZmM37W" ,"ZmMo18W" ,"ZmTx303") ~ genome1,
#                          genome1 %in% c("ZmNC350","ZmNC358") & genome2 %in% c("ZmOh43","ZmOh7b","ZmP39","ZmTx303") ~ genome1,
#                          genome1 == "ZvTIL01" & genome2 == "ZvTIL11" ~ genome1,
#                          .default = genome2
#)
#)
#convergence.shared.stats$genome1_graph<-convergence.shared.stats$genome1_graph %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                                                                     "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                                                     "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                                                     "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,#"ZxTIL25",
#                                                                                                     "ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","TdFL"))
#convergence.shared.stats$genome2_graph<-convergence.shared.stats$genome2_graph %>% factor(levels = c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
#                                                                                                     "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
#                                                                                                     "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
#                                                                                                     "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697"))

#ggplot(data=convergence.shared.stats, aes(x=genome1_graph, y=genome2_graph, fill=mean.allgenes.PercentGeneLost.diff))+
#  facet_wrap(facets=vars(M))+
#  geom_tile()+
#  #  geom_text(aes(label=cnt.NA), size=2, color="darkgrey")+
#  ggtitle("Average percent gene loss difference")+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle=90))+
#  scale_fill_viridis_c(option = "magma", direction = 1, name = "mean")+
#  xlab("")+ylab("")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.avgPercentGeneLossDifference.png",
#       device="png",dpi=300,width=8,height=6)
##Average percent gene loss difference between genomes that are NOT TdFL
#filter(convergence.shared.stats, genome1 != "TdFL" & genome2 != "TdFL") %>% 
#  select(mean.allgenes.PercentGeneLost.diff) %>% pull() %>% mean(na.rm=T)
#Average percent gene loss difference between TdFL and Zea M1
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") & M=="M1") %>% 
#  select(mean.allgenes.PercentGeneLost.diff) %>% pull() %>% mean(na.rm=T)
##Average percent gene loss difference between TdFL and Zea M2
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") & M=="M2") %>% 
#  select(mean.allgenes.PercentGeneLost.diff) %>% pull() %>% mean(na.rm=T)

#Average mean loss sharing absolute between genomes that are NOT TdFL
#filter(convergence.shared.stats, genome1 != "TdFL" & genome2 != "TdFL") %>% 
#  select(mean.LossShareAbs) %>% pull() %>% mean(na.rm=T)
##Average mean loss sharing absolute between TdFL & Zea M1
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") & M=="M1") %>% 
#  select(mean.LossShareAbs ) %>% pull() %>% mean(na.rm=T)
##Average mean loss sharing absolute between TdFL & Zea M2
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") & M=="M2") %>% 
#  select(mean.LossShareAbs) %>% pull() %>% mean(na.rm=T)
##Average mean loss sharing absolute between TdFL & Zea 
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") ) %>% 
#  select(mean.LossShareAbs) %>% pull() %>% mean(na.rm=T)

#Average total loss sharing absolute between genomes that are NOT TdFL
#filter(convergence.shared.stats, genome1 != "TdFL" & genome2 != "TdFL") %>% 
#  select(sum.LossShareAbs) %>% pull() %>% mean(na.rm=T)
##Average total loss sharing absolute between TdFL & Zea M1
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") & M=="M1") %>% 
#  select(sum.LossShareAbs ) %>% pull() %>% mean(na.rm=T)
##Average total loss sharing absolute between TdFL & Zea M2
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") & M=="M2") %>% 
#  select(sum.LossShareAbs) %>% pull() %>% mean(na.rm=T)
##Average total loss sharing absolute between TdFL & Zea 
#filter(convergence.shared.stats, (genome1 == "TdFL" | genome2 == "TdFL") ) %>% 
#  select(sum.LossShareAbs) %>% pull() %>% mean(na.rm=T)

###Better way to visualize???

#ggplot(data=convergence, aes(x=ExonCt, y=PercentGeneLost))+
#  geom_smooth(color="black")+
#  geom_point(aes(color=Genome))+
#  facet_wrap(~M)+
#  theme_bw()+
#  scale_color_manual(values=genome_colors)

#library(ggridges)
#ggplot(convergence.shared,aes(x=PercentGeneLost.diff))+
#  geom_density_ridges(aes(y=genome1,fill=genome1))+
#  theme_minimal()+
#  scale_fill_manual(values = genome_colors)+
#  xlab("Differences")+ylab("")+
#  ggtitle("Percent Gene Lost Differences")

# Are there any instances of a homoeolog being called fractionated but there's NO overlap of which exons were lost? 
#Loss_Pattern_gt1.NoZdZn<-convergence.summary.NoZdZn %>% filter(Number_Unique_Loss_Patterns != 1)
#Loss_Pattern_noOverlap.NoZdZn<-left_join(x=Loss_Pattern_gt1.NoZdZn, y=convergence.shared, by=c("Gene_ID", "M")) %>% filter(Loss_Pattern_Sharing_Per == 0) 
#colnames(Loss_Pattern_noOverlap.NoZdZn)

#ggplot(Loss_Pattern_noOverlap.NoZdZn, aes(x=PercentGeneLost.genome1,y=PercentGeneLost.genome2))+
 # geom_abline(intercept = 0, slope = 1, color = "grey")+
#  geom_point(alpha = 0.4)+
#  theme_bw()+xlab("Percent Gene Lost, genome 1")+ ylab("Percent Gene Lost, genome 2")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/DifferentialFractionation.NoOverlapLossPattern.PercentGeneLost.point.png", device = "png",dpi=300, width = 4, height=4)

#Loss_Pattern_noOverlap.NoZdZn %>% group_by(genome1, genome2, M) %>% count()%>%
# ggplot(aes(x=genome1,y=genome2, fill=n))+
#  geom_tile()+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle = 90))+
#  facet_wrap(~M)+
#  scale_fill_viridis_c(option ="magma", "Counts")+ xlab("")+ylab("")+
#  ggtitle("Total homoeologs showing completely different loss patterns")
#ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/DifferentialFractionation.NoOverlapLossPattern.TotalCounts.heat.png",
#       device = "png", dpi=300, height= 8, width = 11)

#Loss_Pattern_gt1.NoZdZn %>% nrow() #4652

#Loss_Pattern_noOverlap.NoZdZn %>% group_by(genome1, genome2, M) %>% count()%>%ungroup() %>%filter(genome1 != "TdFL" & genome2 != "TdFL") %>% select(n) %>% pull()%>% mean(na.rm=T)
#Loss_Pattern_noOverlap.NoZdZn %>% group_by(genome1, genome2, M) %>% count()%>%ungroup() %>%filter(genome1 == "TdFL" | genome2 == "TdFL") %>% select(n) %>% pull()%>% mean(na.rm=T)
#Loss_Pattern_noOverlap.NoZdZn %>% group_by(genome1, genome2, M) %>% count()%>%ungroup() %>%filter((genome1 == "TdFL" | genome2 == "TdFL") & M== "M1") %>% select(n) %>% pull()%>% mean(na.rm=T)
#Loss_Pattern_noOverlap.NoZdZn %>% group_by(genome1, genome2, M) %>% count()%>%ungroup() %>%filter((genome1 == "TdFL" | genome2 == "TdFL") & M== "M2") %>% select(n) %>% pull()%>% mean(na.rm=T)

#Is the bump at 2 loss patterns strictly differences in fractionation patterns between Tripsacum and Zea? 
#Loss_Pattern_2.NoZdZn<-convergence.summary.NoZdZn %>% filter(Number_Unique_Loss_Patterns == 2)
#left_join(x=Loss_Pattern_2.NoZdZn, y=convergence.shared, by=c("Gene_ID", "M")) %>% filter(Loss_Pattern_Sharing_Abs == 0) %>%
#  group_by(genome1, genome2, M) %>% count() %>%
#  ggplot(aes(x=genome1,y=genome2, fill=n))+
#  geom_tile()+
#  theme_bw()+
#  theme(axis.text.x = element_text(angle = 90))+
#  facet_wrap(~M)+
#  scale_fill_viridis_c("magma")

#left_join(x=Loss_Pattern_2.NoZdZn, y=convergence.shared, by=c("Gene_ID", "M")) %>% filter(Loss_Pattern_Sharing_Abs == 0) %>%
#  group_by(genome1, genome2, M) %>% count() %>% ungroup() %>%filter(genome1 != "TdFL" & genome2 != "TdFL") %>% select(n) %>% pull()%>% mean(na.rm=T)
##average genes shows 2 different loss patterns between Zea genomes
#left_join(x=Loss_Pattern_2.NoZdZn, y=convergence.shared, by=c("Gene_ID", "M")) %>% filter(Loss_Pattern_Sharing_Abs == 0) %>%
#  group_by(genome1, genome2, M) %>% count() %>% ungroup() %>%filter(genome1 == "TdFL" | genome2 == "TdFL") %>% select(n) %>% pull()%>% mean(na.rm=T)
##average genes shows 2 different loss patterns between TdFL and Zea genomes
#left_join(x=Loss_Pattern_2.NoZdZn, y=convergence.shared, by=c("Gene_ID", "M")) %>% filter(Loss_Pattern_Sharing_Abs == 0) %>%
#  group_by(genome1, genome2, M) %>% count() %>% ungroup() %>%filter((genome1 == "TdFL" | genome2 == "TdFL") & M=="M1") %>% select(n) %>% pull()%>% mean(na.rm=T)
##average genes shows 2 different loss patterns between TdFL and Zea genomes M1
#left_join(x=Loss_Pattern_2.NoZdZn, y=convergence.shared, by=c("Gene_ID", "M")) %>% filter(Loss_Pattern_Sharing_Abs == 0) %>%
#  group_by(genome1, genome2, M) %>% count() %>% ungroup() %>%filter((genome1 == "TdFL" | genome2 == "TdFL") & M=="M2") %>% select(n) %>% pull()%>% mean(na.rm=T)
##average genes shows 2 different loss patterns between TdFL and Zea genomes M2
#nrow(Loss_Pattern_2.NoZdZn)

####Revised method for shared.convergence####

#Which genomes have completely different exon loss? Completely the same? Some loss of exons shared but not completely?
#by gene
convergence_shared<-tibble(Gene_ID =NA, Subgenome = NA, 
                           Target_Genome = NA, Query_Genome = NA, Convergence_Category = NA)
for(i in ref.gene){#for each gene
  for(m in c("M1","M2")){ #for each subgenome
    #so long as the gene x subgenome combo isn't in the at least one NA frame
    #i.e. all genomes included have that gene as aligned
    if(nrow(filter(gene_atleastOneNA.NoZnZd, M == m & Gene_ID == i)) == 0){
      #filter the larger convergence data frame to just the gene x subgenome we want to work with
      df<-filter(convergence.NoZnZd, Gene_ID == i & M == m)
      #for each genome as target
      for(g in tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]){
        #pull out the loss pattern for target genome
        pattern1 <- filter(df, Genome == g) %>% select(Loss_Pattern) %>% pull()
        #now for 
        for(q in tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]){
          pattern2 <- filter(df, Genome == q) %>% select(Loss_Pattern) %>% pull()
          if(is_empty(c(pattern1,pattern2))){
            convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                        Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyShared")
          }else{
            #if one of the patterns is empty but not both
            if((is_empty(pattern1) | is_empty(pattern2)) & !is_empty(c(pattern1,pattern2))){
              convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                          Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyDifferent")
              
            }else{
              #if neither pattern is empty
              if(!is_empty(c(pattern1,pattern2))){
                #then do the matching bits
                if(pattern1 == pattern2){ #if the patterns match exactly
                  convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                              Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyShared")
                }else{
                  #figure out which is the longest, but only if both patterns aren't empty
                  if(length(str_split(pattern1, ":",simplify = T)) >= length(str_split(pattern2,":",simplify = T)) & !is_empty(c(pattern1,pattern2))){
                    longest.pattern<-pattern1
                    shortest.pattern<-pattern2
                    #figure out which one is the shortest and longest pattern
                  }else{
                    longest.pattern<-pattern2
                    shortest.pattern<-pattern1
                  }
                  shortest.pattern<-str_split(shortest.pattern,":", simplify=T)
                  #count the number of exons that are detected from shortest pattern in longest pattern
                  counter<-0
                  for(e in 1:length(shortest.pattern)){
                    if(str_detect(longest.pattern, shortest.pattern[1,e])){
                      counter<-counter+1
                    }
                  }
                  #if at least 1 exon was found to be shared between patterns, then somewhat shared
                  if(counter > 0){
                    convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                                Target_Genome = g, Query_Genome = q, Convergence_Category = "SomeShared")
                  }else{ #otherwise it's completely different
                    convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                                Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyDifferent")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
####Need to rerun with the full convergence shared####
#only got through 2304 gene ids (3170 gene x subgenome combinations)
convergence_shared %>% select(Gene_ID, Subgenome) %>% unique() %>% group_by(Subgenome) %>% count()
#2174 for M1 and 995 for M2 (in terms of genes)
convergence_shared<-convergence_shared %>% na.omit()
#how many of the different categories were there:
convergence_shared %>% group_by(Convergence_Category) %>% count()

#heatmap
convergence_shared.stats<-convergence_shared %>% group_by(Subgenome, Target_Genome,Query_Genome,Convergence_Category) %>%
  count()
convergence_shared.stats<-convergence_shared.stats %>% mutate(Percent = case_when(Subgenome == "M1" ~ (n/2174)*100,
                                                                                  Subgenome == "M2" ~ (n/995)*100))
#This would need the fill scales to change I think
ggplot(convergence_shared.stats,aes(x=Target_Genome, y=Query_Genome))+
  geom_tile(aes(fill = Percent))+
  facet_grid(rows = vars(Convergence_Category), cols = vars(Subgenome))+
  xlab("")+ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
  
filter(convergence_shared.stats, Convergence_Category == "CompletelyDifferent") %>%
  ggplot(aes(x=Target_Genome,y=Query_Genome))+
  geom_tile(aes(fill = Percent))+
  facet_wrap(vars(Subgenome))+
  xlab("")+ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Completely Different Exons Fractionated for Same Gene")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.NoZnZd.CDexons.png",
       device = "png", dpi = 300, width = 7, height = 6)
filter(convergence_shared.stats, Convergence_Category == "CompletelyShared") %>%
  ggplot(aes(x=Target_Genome,y=Query_Genome))+
  geom_tile(aes(fill = Percent))+
  facet_wrap(vars(Subgenome))+
  xlab("")+ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Completely Same Exons Fractionated for Same Gene")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.NoZnZd.CSexons.png",
       device = "png", dpi = 300, width = 7, height = 6)
filter(convergence_shared.stats, Convergence_Category == "SomeShared") %>%
  ggplot(aes(x=Target_Genome,y=Query_Genome))+
  geom_tile(aes(fill = Percent))+
  facet_wrap(vars(Subgenome))+
  xlab("")+ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Some Exons Differentially Fractionated for Same Gene")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/convergence.NoZnZd.SSexons.png",
       device = "png", dpi = 300, width = 7, height = 6)

####DIFFERENTIAL RETENTION####
#Here we'll only consider those that are aligned in both M1 and M2
differential_retention<-select(full.fractionation.status, ID, Gene_ID, ends_with("status")) %>% 
  mutate(BothRetainedCnt = rowSums(. == "Both_Retained"),
         M1RetainedCnt = rowSums(. == "M1_Retained"),
         M2RetainedCnt = rowSums(. == "M2_Retained"),
         BothLostCnt = rowSums(. == "Both_Lost")) %>%
  select(contains("ID"),contains("Cnt"))

#Number of exons that are both retained in at least 1 genome and M1/M2 retained in at least 1 genome
filter(differential_retention, BothRetainedCnt > 0 & M1RetainedCnt > 0 & BothLostCnt == 0 & M2RetainedCnt == 0) %>% nrow() #7583
filter(differential_retention, BothRetainedCnt > 0 & M2RetainedCnt > 0 & BothLostCnt == 0 & M1RetainedCnt == 0) %>% nrow() #4867

#Number of Exons that are both lost in at least 1 genome and M1/M2 retained in at least 1 genome
filter(differential_retention, BothLostCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0 & M2RetainedCnt == 0) %>% nrow() #5382
filter(differential_retention, BothLostCnt > 0 & M2RetainedCnt > 0 & BothRetainedCnt == 0 & M1RetainedCnt == 0) %>% nrow() #3539

#Number of Exons that are both kept vs both lost or M1 vs M2 kept in at least 1 genome
filter(differential_retention, BothRetainedCnt > 0 & BothLostCnt > 0 & M2RetainedCnt == 0 & M1RetainedCnt == 0) %>% nrow() #133
filter(differential_retention, BothRetainedCnt == 0 & BothLostCnt == 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0) %>% nrow() #98

#Number of Exons with 3 different fractionation statuses in at least 1 genome
filter(differential_retention, BothRetainedCnt > 0 & BothLostCnt > 0 & M1RetainedCnt > 0 & M2RetainedCnt == 0) %>% nrow() #955
filter(differential_retention, BothRetainedCnt > 0 & BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt == 0) %>% nrow() #486
filter(differential_retention, BothRetainedCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0& BothLostCnt == 0) %>% nrow() #620
filter(differential_retention, BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0) %>% nrow() #621

#Number of Exons with all 4 fractionation statuses in at least 1 genome
filter(differential_retention, BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt > 0) %>% nrow() #615

#Number of Exons where is just 1 fractionation status
filter(differential_retention, BothLostCnt > 0 & M2RetainedCnt == 0 & M1RetainedCnt == 0 & BothRetainedCnt == 0) %>% nrow() #10320
filter(differential_retention, BothLostCnt == 0 & M2RetainedCnt > 0 & M1RetainedCnt == 0 & BothRetainedCnt == 0) %>% nrow() #5966
filter(differential_retention, BothLostCnt == 0 & M2RetainedCnt == 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0) %>% nrow() #15377
filter(differential_retention, BothLostCnt == 0 & M2RetainedCnt == 0 & M1RetainedCnt == 0 & BothRetainedCnt > 0) %>% nrow() #9979

differential_retention.NoZnZd<-select(full.fractionation.status, ID, Gene_ID, ends_with("status")) %>% 
  select(-c(contains("Zn"),"ZdMomo.status","ZdGigi.status")) %>%
  mutate(BothRetainedCnt = rowSums(. == "Both_Retained"),
         M1RetainedCnt = rowSums(. == "M1_Retained"),
         M2RetainedCnt = rowSums(. == "M2_Retained"),
         BothLostCnt = rowSums(. == "Both_Lost")) %>%
  select(contains("ID"),contains("Cnt"))

#Number of exons that are both retained in at least 1 genome and M1/M2 retained in at least 1 genome
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & M1RetainedCnt > 0 & BothLostCnt == 0 & M2RetainedCnt == 0) %>% nrow() #7548
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & M2RetainedCnt > 0 & BothLostCnt == 0 & M1RetainedCnt == 0) %>% nrow() #4769

#Number of Exons that are both lost in at least 1 genome and M1/M2 retained in at least 1 genome
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0 & M2RetainedCnt == 0) %>% nrow() #5246
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M2RetainedCnt > 0 & BothRetainedCnt == 0 & M1RetainedCnt == 0) %>% nrow() #3486

#Number of Exons that are both kept vs both lost or M1 vs M2 kept in at least 1 genome
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & BothLostCnt > 0 & M2RetainedCnt == 0 & M1RetainedCnt == 0) %>% nrow() #133
filter(differential_retention.NoZnZd, BothRetainedCnt == 0 & BothLostCnt == 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0) %>% nrow() #103

#Number of Exons with 3 different fractionation statuses in at least 1 genome
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & BothLostCnt > 0 & M1RetainedCnt > 0 & M2RetainedCnt == 0) %>% nrow() #924
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt == 0) %>% nrow() #448
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothLostCnt == 0) %>% nrow() #620
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0) %>% nrow() #547

#Number of Exons with all 4 fractionation statuses in at least 1 genome
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt > 0) %>% nrow() #555

#Number of Exons where is just 1 fractionation status
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M2RetainedCnt == 0 & M1RetainedCnt == 0 & BothRetainedCnt == 0) %>% nrow() #10360
filter(differential_retention.NoZnZd, BothLostCnt == 0 & M2RetainedCnt > 0 & M1RetainedCnt == 0 & BothRetainedCnt == 0) %>% nrow() #6033
filter(differential_retention.NoZnZd, BothLostCnt == 0 & M2RetainedCnt == 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0) %>% nrow() #15537
filter(differential_retention.NoZnZd, BothLostCnt == 0 & M2RetainedCnt == 0 & M1RetainedCnt == 0 & BothRetainedCnt > 0) %>% nrow() #10178

nrow(differential_retention) #69269

differential_retention.summary<-tibble(NumberOfStatuses = c(1,1,1,1,
                                                            2,2,2,2,2,2,
                                                            3,3,3,3,
                                                            4),
                                       Pattern = c("BothRetained","M1Retained","M2Retained","BothDeleted",
                                                   "BothRetained:M1Retained","BothRetained:M2Retained",
                                                   "BothDeleted:M1Retained","BothDeleted:M2Retained",
                                                   "M1Retained:M2Retained","BothDeleted:BothRetained",
                                                   "BothRetained:M1Retained:M2Retained","BothRetained:M1Retained:BothDeleted",
                                                   "BothRetained:M2Retained:BothDeleted","M1Retained:M2Retained:BothDeleted",
                                                   "BothRetained:M1Retained:M2Retained:BothDeleted"),
                                       Count.AllGenome = c(9979,15377,5966,10320,
                                                           7583,4867,5382,3539,98,133,
                                                           620,955,486,621,
                                                           615),
                                       Count.NoZnZd = c(10178,15537,6033,10360,
                                                        7548,4769,5246,3486,103,133,
                                                        620,924,448,547,
                                                        555))

ggplot(differential_retention.summary, aes(x=NumberOfStatuses, y=Count.AllGenome))+
  geom_point()+
  theme_bw()+
  xlab("Number of Statuses Observed")+ylab("Count of Exons")
  #scale_color_manual(values = c("BothRetained"="#FDE74C","M1Retained"="#7896BF","M2Retained"="#BF4342","BothDeleted"="#2F2F2F",
  #                              "BothRetained:M1Retained"="#80D56F","BothRetained:M2Retained"="#E48220",
  #                              "BothDeleted:M1Retained"="#526683","BothDeleted:M2Retained"="#792929",
  #                              "M1Retained:M2Retained"="#B264EF","BothDeleted:BothRetained"="#A69832",
  #                              "BothRetained:M1Retained:M2Retained"="#5F2F00","BothRetained:M1Retained:BothDeleted"="#548D49",
  #                              "BothRetained:M2Retained:BothDeleted"="#B06418","M1Retained:M2Retained:BothDeleted"="#5D357C",
  #                              "BothRetained:M1Retained:M2Retained:BothDeleted"="#000000"))+
  #theme(legend.text = element_blank())
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/DifferentialFractionation.Counts.png",device="png",dpi=300,height=4,width=4)

differential_retention.summary %>% mutate(Percent.AllGenome = (Count.AllGenome/sum(differential_retention.summary$Count.AllGenome))*100,
                                          Percent.NoZnZd = (Count.NoZnZd/sum(differential_retention.summary$Count.NoZnZd))*100) %>%
  ggplot(aes(x=NumberOfStatuses,y=Percent.AllGenome))+
  geom_point()+
  theme_bw()+
  xlab("Number of Statuses Observed")+ylab("Percentage of Exons")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/DifferentialFractionation.Percentage.png",device="png",dpi=300,height=4,width=4)

#Which genes show all 4 statuses observed?
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt > 0) %>% select(Gene_ID)%>% pull() %>% unique() %>% length()
#235
#Which genes show 3 statuses observed?
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & BothLostCnt > 0 & M1RetainedCnt > 0 & M2RetainedCnt == 0) %>% select(Gene_ID)%>% pull() %>% unique() %>% length() 
#417
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt == 0) %>% select(Gene_ID)%>% pull() %>% unique() %>% length() 
#281
filter(differential_retention.NoZnZd, BothRetainedCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothLostCnt == 0) %>% select(Gene_ID)%>% pull() %>% unique() %>% length() 
#319
filter(differential_retention.NoZnZd, BothLostCnt > 0 & M2RetainedCnt > 0 & M1RetainedCnt > 0 & BothRetainedCnt == 0) %>% select(Gene_ID)%>% pull() %>% unique() %>% length() 
#318


####SPATIAL ORGANIZATION OF FRACTIONATION####
##Does fractionation change by ref chromosome?
filter(ref_Sb313.cds,CHROM == "01")[,4]
long.full.fractionation.status<-long.full.fractionation.status %>% 
  mutate(RefChr = case_when(ID %in% pull(filter(ref_Sb313.cds,CHROM == "01")[,4]) ~ "Chr01",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "02")[,4]) ~ "Chr02",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "03")[,4]) ~ "Chr03",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "04")[,4]) ~ "Chr04",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "05")[,4]) ~ "Chr05",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "06")[,4]) ~ "Chr06",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "07")[,4]) ~ "Chr07",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "08")[,4]) ~ "Chr08",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "09")[,4]) ~ "Chr09",
                            ID %in% pull(filter(ref_Sb313.cds,CHROM == "10")[,4]) ~ "Chr10"))

filter(long.full.fractionation.status, !is.na(RefChr)) %>%
  ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_bw()+
  scale_fill_manual(values = genome_colors)+
  facet_wrap(facets = vars(RefChr))+
  xlab("Number of Exons")+ylab("")+
  ggtitle("Fractionation Status by Exon by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.ByRefChr.AllGenomes.png", device = "png",dpi=500,height=7, width=8)

filter(long.full.fractionation.status, !is.na(RefChr) & !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1")) %>%
  ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_bw()+
  scale_fill_manual(values = genome_colors)+
  facet_wrap(facets = vars(RefChr))+
  xlab("Number of Exons")+ylab("")+
  ggtitle("Fractionation Status by Exon by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.ByRefChr.NoZnZd.png", device = "png",dpi=500,height=7, width=8)


ref_Sb313.cds %>% group_by(CHROM) %>% count()
##For maize meeting, remove the unaligned categories, change to proportion
fractionationStatusByRefChr<-long.full.fractionation.status %>% group_by(Genome,RefChr,Status) %>% count() %>% 
  mutate(Percent = case_when(RefChr == "Chr01" ~ (n/13388)*100,
                             RefChr == "Chr02" ~ (n/8683)*100,
                             RefChr == "Chr03" ~ (n/10580)*100,
                             RefChr == "Chr04" ~ (n/9358)*100,
                             RefChr == "Chr05" ~ (n/2311)*100,
                             RefChr == "Chr06" ~ (n/6169)*100,
                             RefChr == "Chr07" ~ (n/4331)*100,
                             RefChr == "Chr08" ~ (n/2898)*100,
                             RefChr == "Chr09" ~ (n/5743)*100,
                             RefChr == "Chr10" ~ (n/5808)*100))
filter(fractionationStatusByRefChr, !is.na(RefChr)) %>%
  ggplot(aes(x=Percent,y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat="identity")+
  theme_bw()+
  scale_fill_manual(values = genome_colors)+
  facet_wrap(facets = vars(RefChr))+
  xlab("Number of Exons")+ylab("")+
  ggtitle("Fractionation Status by Exon by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.AllCategories.ByRefChr.AllGenomes.png", device = "png",dpi=500,height=7, width=8)


filter(fractionationStatusByRefChr, !is.na(RefChr) & Status %in% c("Both_Lost","M1_Retained","M2_Retained","Both_Retained")) %>%
  ggplot(aes(x=Percent,y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = genome_colors)+
  facet_wrap(facets = vars(RefChr))+
  xlab("Percentage of Exons")+ylab("")+
  ggtitle("Fractionation Status by Exon by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Deleted","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.ByRefChr.AllGenomes.png", device = "png",dpi=500,height=7, width=8)

filter(fractionationStatusByRefChr, !is.na(RefChr) & Status %in% c("Both_Lost","M1_Retained","M2_Retained","Both_Retained") &!Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1")) %>%
  ggplot(aes(x=Percent,y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge(), stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = genome_colors)+
  facet_wrap(facets = vars(RefChr))+
  xlab("Percentage of Exons")+ylab("")+
  ggtitle("Fractionation Status by Exon by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Deleted","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByExon.Percentage.ByRefChr.NoZnZd.png", device = "png",dpi=500,height=7, width=8)

chromosome_colors<-c("Chr01"="#5E0D59", "Chr02"="#8A1465", "Chr03"="#B61B5E", "Chr04"="#E12346", 
                     "Chr05"="#ED524A", "Chr06"="#EE9E7C", "Chr07"="#FFB185", "Chr08"="#FC9E40",
                     "Chr09"="#FEB820", "Chr10"="#FFD53D")
filter(fractionationStatusByRefChr, !is.na(RefChr) & Status %in% c("Both_Retained","Both_Lost","M1_Retained","M2_Retained")) %>%
  ggplot(aes(y=Percent, x=Status, fill=RefChr))+
  geom_violin(position = position_dodge())+
  theme_bw()+
  scale_fill_manual(values = chromosome_colors)+
  ylab("")+xlab("")+ggtitle("Percentage of Exons per Status by Sorghum Ref Chr")+
  scale_x_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Deleted","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 16, face = "bold"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/FractionationStatus.AllGenomes.RefChrColored.Percentage.OnlyAlignedStatuses.png",device="png",dpi=300,
       height = 6, width = 8)
filter(fractionationStatusByRefChr, !is.na(RefChr) & Status %in% c("Both_Retained","Both_Lost","M1_Retained","M2_Retained") & !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1")) %>%
  ggplot(aes(y=Percent, x=Status, fill=RefChr))+
  geom_violin(position = position_dodge())+
  theme_bw()+
  scale_fill_manual(values = chromosome_colors)+
  ylab("")+xlab("")+ggtitle("Percentage of Exons per Status by Sorghum Ref Chr")+
  scale_x_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Deleted","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained"))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size = 16, face = "bold"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/FractionationStatus.NoZnZd.RefChrColored.Percentage.OnlyAlignedStatuses.png",device="png",dpi=300,
       height = 6, width = 8)

filter(fractionationStatusByRefChr, !is.na(RefChr) & Status %in% c("Both_Retained","Both_Lost","M1_Retained","M2_Retained") & !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi")) %>% group_by(RefChr, Status) %>% reframe(mean.percent=mean(Percent, na.rm=T), min.percent=min(Percent, na.rm=T),max.percent=max(Percent, na.rm=T)) %>% view()

fractionationStatusByRefChr$RefChr <-factor(fractionationStatusByRefChr$RefChr)
multcomp::cld(
  multcomp::glht(
    aov(data = filter(fractionationStatusByRefChr, !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1") & Status == "Both_Retained"), 
        formula = Percent ~ RefChr), 
    linfct = multcomp::mcp(RefChr ="Tukey")
    )
  )
#Chr01=A, 02=B, 03=B,04=B,05=C,06=B,07=D,08=E,09=D,10=D
multcomp::cld(
  multcomp::glht(
    aov(data = filter(fractionationStatusByRefChr, !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1") & Status == "M1_Retained"), 
        formula = Percent ~ RefChr), 
    linfct = multcomp::mcp(RefChr ="Tukey")
  )
)
#Chr01=AB, 02=AB, 03=A,04=B,05=C,06=C,07=D,08=E,09=D,10=B

multcomp::cld(
  multcomp::glht(
    aov(data = filter(fractionationStatusByRefChr, !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1") & Status == "M2_Retained"), 
        formula = Percent ~ RefChr), 
    linfct = multcomp::mcp(RefChr ="Tukey")
  )
)
#Chr01=A, 02=A, 03=AB,04=C,05=D,06=C,07=C,08=D,09=D,10=B

multcomp::cld(
  multcomp::glht(
    aov(data = filter(fractionationStatusByRefChr, !Genome %in% c("ZnPI615697","ZdMomo","ZdGigi","ZnPI615697_4to1") & Status == "Both_Lost"), 
        formula = Percent ~ RefChr), 
    linfct = multcomp::mcp(RefChr ="Tukey")
  )
)
#Chr01=A, 02=AB, 03=BC,04=A,05=AB,06=C,07=D,08=BC,09=AB,10=D

sig.diff.byRefChr<-tibble(RefChr_1 = c(rep("Chr01", 36), rep("Chr02", 36), rep("Chr03",36),rep("Chr04",36),rep("Chr05",36),
                                       rep("Chr06", 36), rep("Chr07", 36), rep("Chr08",36),rep("Chr09",36),rep("Chr10",36)),
                          RefChr_2 = c(paste0("Chr0",2:9),"Chr10",paste0("Chr0",2:9),"Chr10",paste0("Chr0",2:9),"Chr10",paste0("Chr0",2:9),"Chr10",
                                       paste0("Chr0",c(1,3:9)),"Chr10",paste0("Chr0",c(1,3:9)),"Chr10",paste0("Chr0",c(1,3:9)),"Chr10",paste0("Chr0",c(1,3:9)),"Chr10",
                                       paste0("Chr0",c(1:2,4:9)),"Chr10",paste0("Chr0",c(1:2,4:9)),"Chr10",paste0("Chr0",c(1:2,4:9)),"Chr10",paste0("Chr0",c(1:2,4:9)),"Chr10",
                                       paste0("Chr0",c(1:3,5:9)),"Chr10",paste0("Chr0",c(1:3,5:9)),"Chr10",paste0("Chr0",c(1:3,5:9)),"Chr10",paste0("Chr0",c(1:3,5:9)),"Chr10",
                                       paste0("Chr0",c(1:4,6:9)),"Chr10",paste0("Chr0",c(1:4,6:9)),"Chr10",paste0("Chr0",c(1:4,6:9)),"Chr10",paste0("Chr0",c(1:4,6:9)),"Chr10",
                                       paste0("Chr0",c(1:5,7:9)),"Chr10",paste0("Chr0",c(1:5,7:9)),"Chr10",paste0("Chr0",c(1:5,7:9)),"Chr10",paste0("Chr0",c(1:5,7:9)),"Chr10",
                                       paste0("Chr0",c(1:6,8:9)),"Chr10",paste0("Chr0",c(1:6,8:9)),"Chr10",paste0("Chr0",c(1:6,8:9)),"Chr10",paste0("Chr0",c(1:6,8:9)),"Chr10",
                                       paste0("Chr0",c(1:7,9)),"Chr10",paste0("Chr0",c(1:7,9)),"Chr10",paste0("Chr0",c(1:7,9)),"Chr10",paste0("Chr0",c(1:7,9)),"Chr10",
                                       paste0("Chr0",c(1:8)),"Chr10",paste0("Chr0",c(1:8)),"Chr10",paste0("Chr0",c(1:8)),"Chr10",paste0("Chr0",c(1:8)),"Chr10",
                                       paste0("Chr0",c(1:9)),paste0("Chr0",c(1:9)),paste0("Chr0",c(1:9)),paste0("Chr0",c(1:9))
                          ),
                          Status = c(rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9),
                                     rep("Both_Retained",9),rep("M1_Retained",9),rep("M2_Retained",9),rep("Both_Deleted",9)
                          ),
                          SignificantlyDifferent=c(rep(T,9),F,F,F,T,T,T,T,T,F,F,F,T,T,T,T,T,T,F,F,T,F,T,F,F,F,T,F,
                                                   T,F,F,T,F,T,T,T,T,F,F,F,T,T,T,T,T,F,F,F,T,T,T,T,T,T,T,F,F,F,F,T,T,F,F,T,
                                                   T,F,F,T,F,T,T,T,T,F,F,T,T,T,T,T,T,T,F,F,T,T,T,T,T,T,F,T,F,T,F,F,T,F,F,T,
                                                   T,F,F,T,F,T,T,T,T,F,F,T,T,T,T,T,T,F,T,T,T,T,F,F,T,T,T,F,F,T,F,T,T,T,F,T,
                                                   T,T,T,T,T,T,T,T,T,T,T,T,T,F,T,T,T,T,T,T,T,T,T,T,F,F,T,F,F,F,F,T,T,F,F,T,
                                                   T,F,F,F,T,T,T,T,T,T,T,T,T,F,T,T,T,T,T,T,T,F,T,F,T,T,T,T,T,F,T,T,T,F,T,T,
                                                   T,T,T,T,T,T,T,F,F,T,T,T,T,T,T,T,F,T,T,T,T,F,T,F,T,T,T,T,T,T,T,T,T,T,T,F,
                                                   T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,T,F,T,T,F,T,T,F,F,T,F,F,T,F,T,
                                                   T,T,T,T,T,T,F,T,F,T,T,T,T,T,T,F,T,T,T,T,T,T,F,T,T,F,T,F,F,F,F,F,T,T,F,T,
                                                   T,T,T,T,T,T,F,T,F,F,F,T,F,T,T,T,T,T,T,T,F,T,T,T,T,T,T,T,T,T,T,T,T,F,T,T)
                          )

ggplot(sig.diff.byRefChr, aes(x=RefChr_1,y=RefChr_2,fill = SignificantlyDifferent))+
  geom_tile()+
  facet_wrap(vars(Status))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))
sig.diff.byRefChr %>% group_by(RefChr_1,RefChr_2,SignificantlyDifferent) %>% count() %>%
  filter(SignificantlyDifferent == T) %>% mutate(percent = (n/4)*100) %>%
ggplot(aes(x=RefChr_1,y=RefChr_2,fill = n))+
  geom_tile()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))

gene_fractionation<-gene_fractionation %>% 
  mutate(RefChr = case_when(Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "01")[,7]) ~ "Chr01",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "02")[,7]) ~ "Chr02",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "03")[,7]) ~ "Chr03",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "04")[,7]) ~ "Chr04",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "05")[,7]) ~ "Chr05",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "06")[,7]) ~ "Chr06",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "07")[,7]) ~ "Chr07",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "08")[,7]) ~ "Chr08",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "09")[,7]) ~ "Chr09",
                            Gene_ID %in% pull(filter(ref_Sb313.cds,CHROM == "10")[,7]) ~ "Chr10"))

long_gene_fractionation<-gene_fractionation %>% select(c(Gene_ID,RefChr,ends_with("status"))) %>% filter(!is.na(RefChr))%>%
  pivot_longer(cols = ends_with("status"), names_to = "status_ID", values_to = "Status") %>%
  filter(Status != "Both_NA") %>%
  mutate(Genome = str_split(status_ID, "[.]", simplify = T)[,1])
long_gene_fractionation$Genome<-long_gene_fractionation$Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                                                                   "ZmHP301" ,"ZmM37W" ,"ZmMo18W" ,"ZmTx303" ,
                                                                                   "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                                                                   "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS"))
long_gene_fractionation$Status<-long_gene_fractionation$Status %>% factor(levels = c("Both_Retained","M1_Retained","M2_Retained", "Both_Lost","M1_Retained:M2_NA","M1_NA:M2_Retained","M1_Lost:M2_NA","M1_NA:M2_Lost","Both_NA"))
ggplot(long_gene_fractionation,aes(y=Status))+
  geom_bar(aes(fill=Genome),position=position_dodge())+
  facet_wrap(facets = vars(RefChr))+
  theme_bw()+
  scale_fill_manual(values=genome_colors)+
  xlab("Number of Genes")+ylab("")+
  ggtitle("Fractionation Status by Gene by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByGene.ByRefChr.AllGenomes.png", dpi=500, device = "png",width = 8,height = 7)

filter(long_gene_fractionation, !Genome %in% c("ZnPI615697","ZdGigi","ZdMomo","ZnPI615697_4to1")) %>% ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome),position=position_dodge())+
  facet_wrap(facets = vars(RefChr))+
  theme_bw()+
  scale_fill_manual(values=genome_colors)+
  xlab("Number of Genes")+ylab("")+
  ggtitle("Fractionation Status by Gene by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/fractionationStatus.ByGene.ByRefChr.NoZnZd.png", dpi=500, device = "png",width = 8,height = 7)

gene.fractionationStatusByRefChr<-long_gene_fractionation %>% group_by(Genome,RefChr,Status) %>% count() %>% 
  mutate(Percent = case_when(RefChr == "Chr01" ~ (n/2343)*100,
                             RefChr == "Chr02" ~ (n/1522)*100,
                             RefChr == "Chr03" ~ (n/1880)*100,
                             RefChr == "Chr04" ~ (n/1545)*100,
                             RefChr == "Chr05" ~ (n/428)*100,
                             RefChr == "Chr06" ~ (n/1107)*100,
                             RefChr == "Chr07" ~ (n/792)*100,
                             RefChr == "Chr08" ~ (n/513)*100,
                             RefChr == "Chr09" ~ (n/1026)*100,
                             RefChr == "Chr10" ~ (n/1001)*100))

gene.fractionationStatusByRefChr %>% filter(Status %in% c("Both_Retained","M1_Retained","M2_Retained","Both_Lost")) %>%
ggplot(aes(x=Percent, y=Status))+
  geom_bar(aes(fill=Genome),position=position_dodge(), stat="identity")+
  facet_wrap(facets = vars(RefChr))+
  theme_bw()+
  scale_fill_manual(values=genome_colors)+
  xlab("Percentage of Genes")+ylab("")+
  ggtitle("Fractionation Status by Gene by Sb Ref Chromosome")+
  scale_y_discrete(labels = c("Both_Retained"="Both Retained","Both_Lost"="Both Lost","M2_Retained"="M2 Retained","M1_Retained"="M1 Retained", "M1_Retained:M2_NA"="M1 Retained, M2 NA","M1_NA:M2_Retained"="M1 NA, M2 Retained","M1_Lost:M2_NA"="M1 Lost, M2 NA","M1_NA:M2_Lost"="M1 NA, M2 Lost","Both_NA"="Both NA"))


#TESTER FILES
#make a bed file of the reference exons for Chr10
#ref_Sb313.cds %>% filter(CHROM == "10") %>% write_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/spatial_fractionation_Chr10/Ref_Exons.WithHeader.bed")
#make a file with the summary of fractionation status for each genome
#filter(full.fractionation.status, ID %in% pull(select(filter(ref_Sb313.cds, CHROM == "10"), ID)))%>% select(c(ID,ends_with("status"))) %>% write_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/spatial_fractionation_Chr10/Fractionation_Status.SummarizedByGenome.tsv")
#make a file with the numeric fractionation calls for each genome
#filter(full.fractionation.status, ID %in% pull(select(filter(ref_Sb313.cds, CHROM == "10"), ID)))%>% select(c(ID,contains(".M"))) %>% write_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/spatial_fractionation_Chr10/Fractation_Status.NumericByM.tsv")

#Plot fractionation amount across genome
##Code from Arun
library(data.table)
library(scales)

#bed<-ref_Sb313.cds %>% filter(CHROM == "10")
#num.Frac<-filter(full.fractionation.status, ID %in% pull(select(filter(ref_Sb313.cds, CHROM == "10"), ID)))%>% select(c(ID,ends_with("status")))
#sum.Frac<-filter(full.fractionation.status, ID %in% pull(select(filter(ref_Sb313.cds, CHROM == "10"), ID)))%>% select(c(ID,contains(".M"))) 
#bed.num.Frac <- merge(bed, num.Frac, by = "ID")
#bed.num.sum.Frac <- merge(bed.num.Frac, sum.Frac, by = "ID")

#myTable <- bed.num.sum.Frac %>%
#  pivot_longer(cols = contains(".M"),
#               names_to = "GenomeM",
#               values_to = "retention") %>%
#  select(CHROM, Start, End, GenomeM, retention)

#myTable$retention <- factor(myTable$retention, levels = c("0", "1"))
#ggplot(
#  myTable) +
#  geom_segment(aes(
#    x = Start,
#    xend = End,
#    y = GenomeM,
#    yend = GenomeM,
#    color = retention),
#    size = 1.0
#  ) +
#  labs(x = "Genomic Position", y = "Chromosome") +
#  ggtitle("Exon Regions") +
#  expand_limits(x = c(0, NA), y = c(0, NA)) +
#  scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
#  theme_minimal()

#myTable$retention <- as.numeric(myTable$retention)
#myTable$retention <- myTable$retention - 1
#newTable <- myTable %>% 
#  mutate(M = gsub(".*\\.M", "M", GenomeM)) %>% 
#  select(CHROM, Start, End, M, retention) %>%
#  group_by(CHROM, Start, End, M) %>%
#  summarise(score = sum(retention, na.rm = TRUE))

#ggplot(newTable) +
#  geom_segment(aes(
#    x = Start,
#    xend = End,
#    y = 0,
#   yend = score,
#    color = M),
#    size = 1.0
#  ) +
#  labs(x = "Genomic Position", y = "Score") +
#  ggtitle("frequency of exon deletion") +
#  expand_limits(x = c(0, NA), y = c(0, NA)) +
#  scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
#  theme_bw() + facet_wrap(~M, ncol = 1)+theme(legend.position = "none")

#Let's turn it into a function
exonDelFreq_byRefChr<-function(chr){ #chr="10"
  bed<-ref_Sb313.cds %>% filter(CHROM == chr)
  num.Frac<-filter(full.fractionation.status, ID %in% pull(select(filter(ref_Sb313.cds, CHROM == chr), ID)))%>% select(c(ID,ends_with("status")))
  sum.Frac<-filter(full.fractionation.status, ID %in% pull(select(filter(ref_Sb313.cds, CHROM == chr), ID)))%>% select(c(ID,contains(".M"))) 
  bed.num.Frac <- merge(bed, num.Frac, by = "ID")
  bed.num.sum.Frac <- merge(bed.num.Frac, sum.Frac, by = "ID") 
  myTable <- bed.num.sum.Frac %>%
    pivot_longer(cols = contains(".M"),
                 names_to = "GenomeM",
                 values_to = "retention") %>%
    select(CHROM, Start, End, GenomeM, retention)
  myTable$retention <- factor(myTable$retention, levels = c("0", "1"))
  myTable$retention <- as.numeric(myTable$retention)
  myTable$retention <- myTable$retention - 1
  newTable <- myTable %>% 
    mutate(M = gsub(".*\\.M", "M", GenomeM)) %>% 
    select(CHROM, Start, End, M, retention) %>%
    group_by(CHROM, Start, End, M) %>%
    summarise(score = sum(retention, na.rm = TRUE))
  
  p<-ggplot(newTable) +
    geom_segment(aes(
      x = Start,
      xend = End,
      y = 0,
      yend = score,
      color = M),
      size = 1.0
    ) +
    labs(x = "Genomic Position", y = "Score") +
    ggtitle(paste("Frequency of Exon Deletion Across Sb Chr",chr)) +
    expand_limits(x = c(0, NA), y = c(0, NA)) +
    scale_x_continuous(labels = unit_format(unit = "Mb", scale = 1e-6)) +
    scale_color_manual(values=subgenome_colors)+
    theme_bw() + facet_wrap(~M, ncol = 1)+theme(legend.position = "none")
  return(p)
}

exonDelFreq_byRefChr("10")
exonDelFreq_byRefChr("01")

for(i in c("01","02","03","04","05","06","07","08","09","10")){
  exonDelFreq_byRefChr(i)
  ggsave(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/freqExonDel.ByPhysPos.Sb.",i,".png"),
         device="png",dpi=300, width=6, height=5)
}

#Test if there's differences in retention status by Ref_Chr and Genome
#Chi-squared test
long_gene_fractionation %>% ungroup() %>% group_by(RefChr, Status) %>% count()
long_gene_fractionation %>% ungroup() %>% group_by(RefChr, Status) %>% count() %>% pivot_wider(names_from = RefChr, values_from = n, id_cols = Status)
table(long_gene_fractionation$RefChr, long_gene_fractionation$Status)

chisq.test(long_gene_fractionation$RefChr, long_gene_fractionation$Status)
#data:  long_gene_fractionation$RefChr and long_gene_fractionation$Status
#X-squared = 18700, df = 63, p-value < 2.2e-16

#to do the post-hoc test
devtools::install_github("ebbertd/chisq.posthoc.test")

#chisq.posthoc.test::chisq.posthoc.test(table(long_gene_fractionation$RefChr, long_gene_fractionation$Status))
#Not working

#TEST FOR RETAINED#
filter(long.full.fractionation.status, Status == "Both_Retained") %>% .[,c("Genome","RefChr")] %>% table()
filter(long.full.fractionation.status, Status == "Both_Retained") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.test() 
#X-squared = 4181.9, df = 342, p-value < 2.2e-16
#filter(long.full.fractionation.status, Status == "Both_Retained") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.posthoc.test::chisq.posthoc.test()%>%
#  filter(Value == "p values")


filter(long.full.fractionation.status, Status == "Both_Retained" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table() %>%chisq.test() 
#X-squared = 2710.6, df = 315, p-value < 2.2e-16
filter(long.full.fractionation.status, Status == "Both_Retained" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table()%>%chisq.posthoc.test::chisq.posthoc.test()%>%
  filter(Value == "p values")
#Retention was significant for

#TEST FOR BOTH LOST#
filter(long.full.fractionation.status, Status == "Both_Lost") %>% .[,c("Genome","RefChr")] %>% table()
filter(long.full.fractionation.status, Status == "Both_Lost") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.test() 
#X-squared = 3333.6, df = 342, p-value < 2.2e-16
#filter(long.full.fractionation.status, Status == "Both_Lost") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.posthoc.test::chisq.posthoc.test()%>%
#  filter(Value == "p values")
#Retention was significant for

filter(long.full.fractionation.status, Status == "Both_Lost" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table() %>%chisq.test() 
#X-squared = 2041.4, df = 315, p-value < 2.2e-16
filter(long.full.fractionation.status, Status == "Both_Lost" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table()%>%chisq.posthoc.test::chisq.posthoc.test()%>%
  filter(Value == "p values")
#Retention was significant for

#TEST FOR M1 RETAINED#
filter(long.full.fractionation.status, Status == "M1_Retained") %>% .[,c("Genome","RefChr")] %>% table()
filter(long.full.fractionation.status, Status == "M1_Retained") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.test() 
#X-squared = 4581.5, df = 342, p-value < 2.2e-16
#filter(long.full.fractionation.status, Status == "M1_Retained") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.posthoc.test::chisq.posthoc.test()%>%
#  filter(Value == "p values")
#Retention was significant for

filter(long.full.fractionation.status, Status == "M1_Retained" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table() %>%chisq.test() 
#X-squared =2831.2, df = 315, p-value < 2.2e-16
#filter(long.full.fractionation.status, Status == "M1_Retained" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table()%>%chisq.posthoc.test::chisq.posthoc.test()%>%
#  filter(Value == "p values")
#Retention was significant for

#TEST FOR M2 RETAINED
filter(long.full.fractionation.status, Status == "M2_Retained") %>% .[,c("Genome","RefChr")] %>% table()
filter(long.full.fractionation.status, Status == "M2_Retained") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.test()
#X-squared = 2526.6, df = 342, p-value < 2.2e-16
#filter(long.full.fractionation.status, Status == "M2_Retained") %>% .[,c("Genome","RefChr")] %>% table() %>%chisq.posthoc.test::chisq.posthoc.test()%>%
#  filter(Value == "p values")
#Retention was significant for

filter(long.full.fractionation.status, Status == "M2_Retained" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table() %>%chisq.test() 
# X-squared = 2526.6, df = 342, p-value < 2.2e-16
#filter(long.full.fractionation.status, Status == "M2_Retained" & !Genome %in% c("ZdMomo","ZdGigi","ZnPI615697")) %>% .[,c("Genome","RefChr")] %>% droplevels()%>% table()%>%chisq.posthoc.test::chisq.posthoc.test()%>%
#  filter(Value == "p values")
#Retention was significant for

filter(long.full.fractionation.status) %>% group_by(RefChr,Status) %>% count()

####ESTIMATE TIMING OF FRACTIONATION####
#For a deletion
#if it is shared by all genomes of M1 | M2 --> basal deletion
#if it is private to one genome of M1 | M2 --> species deletion
#if it is only in Zea genomes of M1 | M2 --> genus deletion
#if it is only in Zv and Zm of M1 | M2 --> sister deletion
#if it is any other sharing of M1 | M2 --> paraphyly deletion

timing_categories<-read_csv("/work/LAS/mhufford-lab/snodgras/Fractionation/Timing_Categories.csv")
timing_categories[nrow(timing_categories) + 1, 3:length(timing_categories)]<-0
timing_categories[nrow(timing_categories), "Age"]<-"CompletelyRetained"
timing_categories[nrow(timing_categories),"Node"]<-"N0_Retained"
timing_categories_withoutZdZn<-filter(timing_categories, !Node %in% c("N2_Zn","N2_ZdZm")) %>% select(-c(ZnPI615697, ZdGigi,ZdMomo,ZnPI615697_4to1))

full.fractionation.status$M1.Timing<-NA
full.fractionation.status$M1.TimingNode<-NA
full.fractionation.status$M2.Timing<-NA
full.fractionation.status$M2.TimingNode<-NA

for(i in 1:nrow(full.fractionation.status)){
  if(all(!is.na(select(full.fractionation.status, ends_with(".M1"))[i,]))){ #If all is not NA for M1
    df<-filter(timing_categories, TdFL == pull(full.fractionation.status[i,"TdFL.M1"]) &
                 TdKS == pull(full.fractionation.status[i,"TdKS.M1"]) &
                 ZdGigi_4to1 == pull(full.fractionation.status[i,paste0("ZdGigi_4to1",".M1")]) &
                 ZdMomo_4to1 == pull(full.fractionation.status[i,paste0("ZdMomo_4to1",".M1")]) &
                 ZdGigi == pull(full.fractionation.status[i,paste0("ZdGigi",".M1")]) &
                 ZdMomo == pull(full.fractionation.status[i,paste0("ZdMomo",".M1")]) &
                 ZhRIMHU001 == pull(full.fractionation.status[i,paste0("ZhRIMHU001",".M1")]) &
                 ZmB73 == pull(full.fractionation.status[i,paste0("ZmB73",".M1")]) &
                 ZmB97 == pull(full.fractionation.status[i,paste0("ZmB97",".M1")]) &
                 ZmCML103 == pull(full.fractionation.status[i,paste0("ZmCML103",".M1")]) &
                 ZmCML228 == pull(full.fractionation.status[i,paste0("ZmCML228",".M1")]) &
                 ZmCML247 == pull(full.fractionation.status[i,paste0("ZmCML247",".M1")]) &
                 ZmCML277 == pull(full.fractionation.status[i,paste0("ZmCML277",".M1")]) &
                 ZmCML322 == pull(full.fractionation.status[i,paste0("ZmCML322",".M1")]) &
                 ZmCML333 == pull(full.fractionation.status[i,paste0("ZmCML333",".M1")]) &
                 ZmCML52 == pull(full.fractionation.status[i,paste0("ZmCML52",".M1")]) &
                 ZmCML69 == pull(full.fractionation.status[i,paste0("ZmCML69",".M1")]) &
                 ZmHP301 == pull(full.fractionation.status[i,paste0("ZmHP301",".M1")]) &
                 ZmIL14H == pull(full.fractionation.status[i,paste0("ZmIL14H",".M1")]) &
                 ZmKi11 == pull(full.fractionation.status[i,paste0("ZmKi11",".M1")]) &
                 ZmKi3 == pull(full.fractionation.status[i,paste0("ZmKi3",".M1")]) &
                 ZmKy21 == pull(full.fractionation.status[i,paste0("ZmKy21",".M1")]) &
                 ZmM162W == pull(full.fractionation.status[i,paste0("ZmM162W",".M1")]) &
                 
                 ZmMS71 == pull(full.fractionation.status[i,paste0("ZmMS71",".M1")]) &
                 ZmNC350 == pull(full.fractionation.status[i,paste0("ZmNC350",".M1")]) &
                 ZmNC358 == pull(full.fractionation.status[i,paste0("ZmNC358",".M1")]) &
                 ZmOh43 == pull(full.fractionation.status[i,paste0("ZmOh43",".M1")]) &
                 ZmOh7b == pull(full.fractionation.status[i,paste0("ZmOh7b",".M1")]) &
                 ZmP39 == pull(full.fractionation.status[i,paste0("ZmP39",".M1")]) &
                 
                 ZmTzi8 == pull(full.fractionation.status[i,paste0("ZmTzi8",".M1")]) &
                 ZnPI615697 == pull(full.fractionation.status[i,paste0("ZnPI615697",".M1")]) &
                 ZnPI615697_4to1 == pull(full.fractionation.status[i,paste0("ZnPI615697_4to1",".M1")]) &
                 ZvTIL01 == pull(full.fractionation.status[i,paste0("ZvTIL01",".M1")]) &
                 ZvTIL11 == pull(full.fractionation.status[i,paste0("ZvTIL11",".M1")]) &
                 ZxTIL18 == pull(full.fractionation.status[i,paste0("ZxTIL18",".M1")]) &
                 ZxTIL25 == pull(full.fractionation.status[i,paste0("ZxTIL25",".M1")]) )
    if(nrow(df) > 0){ #if there are no NAs and there's a timing pattern match:
      full.fractionation.status$M1.Timing[i]<-df %>% select("Age") %>% pull()
      full.fractionation.status$M1.TimingNode[i]<-df %>% select("Node") %>% pull()
    }else{#else if there's no match to a timing pattern:
      full.fractionation.status$M1.Timing[i]<-"paraphyly"
      full.fractionation.status$M1.TimingNode[i]<-"paraphyly"
    }
  } 
  if(all(is.na(select(full.fractionation.status, ends_with(".M1"))[i,]))){ #if everything is NA
    full.fractionation.status$M1.Timing[i]<-"Completely_Unaligned"
    full.fractionation.status$M1.TimingNode[i]<-NA
  }
  if(any(is.na(select(full.fractionation.status, ends_with(".M1"))[i,])) & !all(is.na(select(full.fractionation.status, ends_with(".M1"))[i,]))){ #if anything is NA, but not everything
    full.fractionation.status$M1.Timing[i]<-"Some_Unaligned"
    full.fractionation.status$M1.TimingNode[i]<-NA
  }
  #Same loop now with M2
  if(all(!is.na(select(full.fractionation.status, ends_with(".M2"))[i,]))){
    df<-filter(timing_categories, TdFL == pull(full.fractionation.status[i,"TdFL.M2"]) &
                 TdKS == pull(full.fractionation.status[i,"TdKS.M2"]) &
                 ZdGigi_4to1 == pull(full.fractionation.status[i,paste0("ZdGigi_4to1",".M2")]) &
                 ZdMomo_4to1 == pull(full.fractionation.status[i,paste0("ZdMomo_4to1",".M2")]) &
                 ZdGigi == pull(full.fractionation.status[i,paste0("ZdGigi",".M2")]) &
                 ZdMomo == pull(full.fractionation.status[i,paste0("ZdMomo",".M2")]) &
                 ZhRIMHU001 == pull(full.fractionation.status[i,paste0("ZhRIMHU001",".M2")]) &
                 ZmB73 == pull(full.fractionation.status[i,paste0("ZmB73",".M2")]) &
                 ZmB97 == pull(full.fractionation.status[i,paste0("ZmB97",".M2")]) &
                 ZmCML103 == pull(full.fractionation.status[i,paste0("ZmCML103",".M2")]) &
                 ZmCML228 == pull(full.fractionation.status[i,paste0("ZmCML228",".M2")]) &
                 ZmCML247 == pull(full.fractionation.status[i,paste0("ZmCML247",".M2")]) &
                 ZmCML277 == pull(full.fractionation.status[i,paste0("ZmCML277",".M2")]) &
                 ZmCML322 == pull(full.fractionation.status[i,paste0("ZmCML322",".M2")]) &
                 ZmCML333 == pull(full.fractionation.status[i,paste0("ZmCML333",".M2")]) &
                 ZmCML52 == pull(full.fractionation.status[i,paste0("ZmCML52",".M2")]) &
                 ZmCML69 == pull(full.fractionation.status[i,paste0("ZmCML69",".M2")]) &
                 ZmHP301 == pull(full.fractionation.status[i,paste0("ZmHP301",".M2")]) &
                 ZmIL14H == pull(full.fractionation.status[i,paste0("ZmIL14H",".M2")]) &
                 ZmKi11 == pull(full.fractionation.status[i,paste0("ZmKi11",".M2")]) &
                 ZmKi3 == pull(full.fractionation.status[i,paste0("ZmKi3",".M2")]) &
                 ZmKy21 == pull(full.fractionation.status[i,paste0("ZmKy21",".M2")]) &
                 ZmM162W == pull(full.fractionation.status[i,paste0("ZmM162W",".M2")]) &
                 
                 ZmMS71 == pull(full.fractionation.status[i,paste0("ZmMS71",".M2")]) &
                 ZmNC350 == pull(full.fractionation.status[i,paste0("ZmNC350",".M2")]) &
                 ZmNC358 == pull(full.fractionation.status[i,paste0("ZmNC358",".M2")]) &
                 ZmOh43 == pull(full.fractionation.status[i,paste0("ZmOh43",".M2")]) &
                 ZmOh7b == pull(full.fractionation.status[i,paste0("ZmOh7b",".M2")]) &
                 ZmP39 == pull(full.fractionation.status[i,paste0("ZmP39",".M2")]) &
                 
                 ZmTzi8 == pull(full.fractionation.status[i,paste0("ZmTzi8",".M2")]) &
                 ZnPI615697 == pull(full.fractionation.status[i,paste0("ZnPI615697",".M2")]) &
                 ZnPI615697_4to1 == pull(full.fractionation.status[i,paste0("ZnPI615697_4to1",".M2")]) &
                 ZvTIL01 == pull(full.fractionation.status[i,paste0("ZvTIL01",".M2")]) &
                 ZvTIL11 == pull(full.fractionation.status[i,paste0("ZvTIL11",".M2")]) &
                 ZxTIL18 == pull(full.fractionation.status[i,paste0("ZxTIL18",".M2")]) &
                 ZxTIL25 == pull(full.fractionation.status[i,paste0("ZxTIL25",".M2")]) )
    if(nrow(df) > 0){
      full.fractionation.status$M2.Timing[i]<-df %>% select("Age") %>% pull()
      full.fractionation.status$M2.TimingNode[i]<-df %>% select("Node") %>% pull()
    }else{
      full.fractionation.status$M2.Timing[i]<-"paraphyly"
      full.fractionation.status$M2.TimingNode[i]<-"paraphyly"
    }
  }
  if(all(is.na(select(full.fractionation.status, ends_with(".M2"))[i,]))){ #if everything is NA
    full.fractionation.status$M2.Timing[i]<-"Completely_Unaligned"
    full.fractionation.status$M2.TimingNode[i]<-NA
  }
  if(any(is.na(select(full.fractionation.status, ends_with(".M2"))[i,])) & !all(is.na(select(full.fractionation.status, ends_with(".M2"))[i,]))){ #if anything is NA, but not everything
    full.fractionation.status$M2.Timing[i]<-"Some_Unaligned"
    full.fractionation.status$M2.TimingNode[i]<-NA
  }
}

full.fractionation.status$M1.Timing.NoZnZd<-NA
full.fractionation.status$M1.TimingNode.NoZnZd<-NA
full.fractionation.status$M2.Timing.NoZnZd<-NA
full.fractionation.status$M2.TimingNode.NoZnZd<-NA

for(i in 1:nrow(full.fractionation.status)){
  if(all(!is.na(select(full.fractionation.status, ends_with(".M1"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,]))){
    df<-filter(timing_categories_withoutZdZn, TdFL == pull(full.fractionation.status[i,"TdFL.M1"]) &
                 TdKS == pull(full.fractionation.status[i,"TdKS.M1"]) &
                 ZdGigi_4to1 == pull(full.fractionation.status[i,paste0("ZdGigi_4to1",".M1")]) &
                 ZdMomo_4to1 == pull(full.fractionation.status[i,paste0("ZdMomo_4to1",".M1")]) &
                 ZhRIMHU001 == pull(full.fractionation.status[i,paste0("ZhRIMHU001",".M1")]) &
                 ZmB73 == pull(full.fractionation.status[i,paste0("ZmB73",".M1")]) &
                 ZmB97 == pull(full.fractionation.status[i,paste0("ZmB97",".M1")]) &
                 ZmCML103 == pull(full.fractionation.status[i,paste0("ZmCML103",".M1")]) &
                 ZmCML228 == pull(full.fractionation.status[i,paste0("ZmCML228",".M1")]) &
                 ZmCML247 == pull(full.fractionation.status[i,paste0("ZmCML247",".M1")]) &
                 ZmCML277 == pull(full.fractionation.status[i,paste0("ZmCML277",".M1")]) &
                 ZmCML322 == pull(full.fractionation.status[i,paste0("ZmCML322",".M1")]) &
                 ZmCML333 == pull(full.fractionation.status[i,paste0("ZmCML333",".M1")]) &
                 ZmCML52 == pull(full.fractionation.status[i,paste0("ZmCML52",".M1")]) &
                 ZmCML69 == pull(full.fractionation.status[i,paste0("ZmCML69",".M1")]) &
                 ZmHP301 == pull(full.fractionation.status[i,paste0("ZmHP301",".M1")]) &
                 ZmIL14H == pull(full.fractionation.status[i,paste0("ZmIL14H",".M1")]) &
                 ZmKi11 == pull(full.fractionation.status[i,paste0("ZmKi11",".M1")]) &
                 ZmKi3 == pull(full.fractionation.status[i,paste0("ZmKi3",".M1")]) &
                 ZmKy21 == pull(full.fractionation.status[i,paste0("ZmKy21",".M1")]) &
                 ZmM162W == pull(full.fractionation.status[i,paste0("ZmM162W",".M1")]) &
                 
                 
                 ZmMS71 == pull(full.fractionation.status[i,paste0("ZmMS71",".M1")]) &
                 ZmNC350 == pull(full.fractionation.status[i,paste0("ZmNC350",".M1")]) &
                 ZmNC358 == pull(full.fractionation.status[i,paste0("ZmNC358",".M1")]) &
                 ZmOh43 == pull(full.fractionation.status[i,paste0("ZmOh43",".M1")]) &
                 ZmOh7b == pull(full.fractionation.status[i,paste0("ZmOh7b",".M1")]) &
                 ZmP39 == pull(full.fractionation.status[i,paste0("ZmP39",".M1")]) &
                 
                 ZmTzi8 == pull(full.fractionation.status[i,paste0("ZmTzi8",".M1")]) &
                 #ZnPI615697 == pull(full.fractionation.status[i,paste0("ZnPI615697",".M1")]) &
                 ZvTIL01 == pull(full.fractionation.status[i,paste0("ZvTIL01",".M1")]) &
                 ZvTIL11 == pull(full.fractionation.status[i,paste0("ZvTIL11",".M1")]) &
                 ZxTIL18 == pull(full.fractionation.status[i,paste0("ZxTIL18",".M1")]) &
                 ZxTIL25 == pull(full.fractionation.status[i,paste0("ZxTIL25",".M1")]) )
    if(nrow(df) > 0){
      full.fractionation.status$M1.Timing.NoZnZd[i]<-df %>% select("Age") %>% pull()
      full.fractionation.status$M1.TimingNode.NoZnZd[i]<-df %>% select("Node") %>% pull()
    }else{
      full.fractionation.status$M1.Timing.NoZnZd[i]<-"paraphyly"
      full.fractionation.status$M1.TimingNode.NoZnZd[i]<-"paraphyly"
    }
  } 
  if(all(is.na(select(full.fractionation.status, ends_with(".M1"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,]))){ #if everything is NA
    full.fractionation.status$M1.Timing.NoZnZd[i]<-"Completely_Unaligned"
    full.fractionation.status$M1.TimingNode.NoZnZd[i]<-NA
  }
  if(any(is.na(select(full.fractionation.status, ends_with(".M1"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,])) & !all(is.na(select(full.fractionation.status, ends_with(".M1"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,]))){ #if anything is NA, but not everything
    full.fractionation.status$M1.Timing.NoZnZd[i]<-"Some_Unaligned"
    full.fractionation.status$M1.TimingNode.NoZnZd[i]<-NA
  }
  #Same loop now with M2
  if(all(!is.na(select(full.fractionation.status, ends_with(".M2"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,]))){
    df<-filter(timing_categories_withoutZdZn, TdFL == pull(full.fractionation.status[i,"TdFL.M2"]) &
                 TdKS == pull(full.fractionation.status[i,"TdKS.M2"]) &
                 ZdGigi_4to1 == pull(full.fractionation.status[i,paste0("ZdGigi_4to1",".M2")]) &
                 ZdMomo_4to1 == pull(full.fractionation.status[i,paste0("ZdMomo_4to1",".M2")]) &
                 ZhRIMHU001 == pull(full.fractionation.status[i,paste0("ZhRIMHU001",".M2")]) &
                 ZmB73 == pull(full.fractionation.status[i,paste0("ZmB73",".M2")]) &
                 ZmB97 == pull(full.fractionation.status[i,paste0("ZmB97",".M2")]) &
                 ZmCML103 == pull(full.fractionation.status[i,paste0("ZmCML103",".M2")]) &
                 ZmCML228 == pull(full.fractionation.status[i,paste0("ZmCML228",".M2")]) &
                 ZmCML247 == pull(full.fractionation.status[i,paste0("ZmCML247",".M2")]) &
                 ZmCML277 == pull(full.fractionation.status[i,paste0("ZmCML277",".M2")]) &
                 ZmCML322 == pull(full.fractionation.status[i,paste0("ZmCML322",".M2")]) &
                 ZmCML333 == pull(full.fractionation.status[i,paste0("ZmCML333",".M2")]) &
                 ZmCML52 == pull(full.fractionation.status[i,paste0("ZmCML52",".M2")]) &
                 ZmCML69 == pull(full.fractionation.status[i,paste0("ZmCML69",".M2")]) &
                 ZmHP301 == pull(full.fractionation.status[i,paste0("ZmHP301",".M2")]) &
                 ZmIL14H == pull(full.fractionation.status[i,paste0("ZmIL14H",".M2")]) &
                 ZmKi11 == pull(full.fractionation.status[i,paste0("ZmKi11",".M2")]) &
                 ZmKi3 == pull(full.fractionation.status[i,paste0("ZmKi3",".M2")]) &
                 ZmKy21 == pull(full.fractionation.status[i,paste0("ZmKy21",".M2")]) &
                 ZmM162W == pull(full.fractionation.status[i,paste0("ZmM162W",".M2")]) &
                 
                 
                 ZmMS71 == pull(full.fractionation.status[i,paste0("ZmMS71",".M2")]) &
                 ZmNC350 == pull(full.fractionation.status[i,paste0("ZmNC350",".M2")]) &
                 ZmNC358 == pull(full.fractionation.status[i,paste0("ZmNC358",".M2")]) &
                 ZmOh43 == pull(full.fractionation.status[i,paste0("ZmOh43",".M2")]) &
                 ZmOh7b == pull(full.fractionation.status[i,paste0("ZmOh7b",".M2")]) &
                 ZmP39 == pull(full.fractionation.status[i,paste0("ZmP39",".M2")]) &
                 
                 ZmTzi8 == pull(full.fractionation.status[i,paste0("ZmTzi8",".M2")]) &
                 #ZnPI615697 == pull(full.fractionation.status[i,paste0("ZnPI615697",".M2")]) &
                 ZvTIL01 == pull(full.fractionation.status[i,paste0("ZvTIL01",".M2")]) &
                 ZvTIL11 == pull(full.fractionation.status[i,paste0("ZvTIL11",".M2")]) &
                 ZxTIL18 == pull(full.fractionation.status[i,paste0("ZxTIL18",".M2")]) &
                 ZxTIL25 == pull(full.fractionation.status[i,paste0("ZxTIL25",".M2")]) )
    if(nrow(df) > 0){
      full.fractionation.status$M2.Timing.NoZnZd[i]<-df %>% select("Age") %>% pull()
      full.fractionation.status$M2.TimingNode.NoZnZd[i]<-df %>% select("Node") %>% pull()
    }else{
      full.fractionation.status$M2.Timing.NoZnZd[i]<-"paraphyly"
      full.fractionation.status$M2.TimingNode.NoZnZd[i]<-"paraphyly"
    }
  }
  if(all(is.na(select(full.fractionation.status, ends_with(".M2"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,]))){
    full.fractionation.status$M2.Timing.NoZnZd[i]<-"Completely_Unaligned"
    full.fractionation.status$M2.TimingNode.NoZnZd[i]<-NA
  }
  if(any(is.na(select(full.fractionation.status, ends_with(".M2"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,])) & !all(is.na(select(full.fractionation.status, ends_with(".M2"),-c(contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))[i,]))){ #if anything is NA, but not everything
    full.fractionation.status$M2.Timing.NoZnZd[i]<-"Some_Unaligned"
    full.fractionation.status$M2.TimingNode.NoZnZd[i]<-NA
  }
}

#Make timing a factor 
full.fractionation.status$M1.Timing <- full.fractionation.status$M1.Timing %>% factor(c("basal","genus","genus_sister","species","subspecies_sister1","subspecies_sister2","subspecies_private","accession_sister","private","CompletelyRetained","paraphyly","Some_Unaligned","Completely_Unaligned"))
full.fractionation.status$M2.Timing <- full.fractionation.status$M2.Timing %>% factor(c("basal","genus","genus_sister","species","subspecies_sister1","subspecies_sister2","subspecies_private","accession_sister","private","CompletelyRetained","paraphyly","Some_Unaligned","Completely_Unaligned"))
full.fractionation.status$M1.Timing.NoZnZd <- full.fractionation.status$M1.Timing.NoZnZd %>% factor(c("basal","genus","genus_sister","species","subspecies_sister1","subspecies_sister2","subspecies_private","accession_sister","private","CompletelyRetained","paraphyly","Some_Unaligned","Completely_Unaligned"))
full.fractionation.status$M2.Timing.NoZnZd <- full.fractionation.status$M2.Timing.NoZnZd %>% factor(c("basal","genus","genus_sister","species","subspecies_sister1","subspecies_sister2","subspecies_private","accession_sister","private","CompletelyRetained","paraphyly","Some_Unaligned","Completely_Unaligned"))

tibble(Timing = c(full.fractionation.status$M1.Timing, full.fractionation.status$M2.Timing),
       M = c(rep("M1",nrow(full.fractionation.status)),rep("M2",nrow(full.fractionation.status)))) %>%
  filter(!is.na(Timing))%>%
  ggplot(aes(y=Timing))+
  geom_bar(aes(fill = M),position=position_dodge())+
  scale_fill_manual(values=subgenome_colors)+
  theme_minimal()+xlab("Count by Exon")+
  ggtitle("Timing using Exon, All Genomes")+
  scale_y_discrete(labels = c("basal" = "Basal", "genus"="Genus","genus_sister"="Genus, sister","species"="Species","subspecies_private"="Subspecies, private","accession_sister"="Accession, sister","private"="Private","CompletelyRetained"="Completely Retained","paraphyly"="Paraphyly", "Some_Unaligned" = "Some Unaligned","Completely_Unaligned"="Completely Unaligned"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.generalcategories.allgenomes.png", device="png",dpi=300, width = 3, height = 3)

tibble(Timing = c(full.fractionation.status$M1.Timing, full.fractionation.status$M2.Timing),
       M = c(rep("M1",nrow(full.fractionation.status)),rep("M2",nrow(full.fractionation.status)))) %>%
  filter(!is.na(Timing) & !Timing %in% c("Some_Unaligned","Completely_Unaligned"))%>%
  ggplot(aes(y=Timing))+
  geom_bar(aes(fill = M),position=position_dodge())+
  scale_fill_manual(values=subgenome_colors)+
  theme_minimal()+xlab("Count by Exon")+
  ggtitle("Timing using Exon, All Genomes")+
  scale_y_discrete(labels = c("basal" = "Basal", "genus"="Genus","genus_sister"="Genus, sister","species"="Species","subspecies_private"="Subspecies, private","accession_sister"="Accession, sister","private"="Private","CompletelyRetained"="Completely Retained","paraphyly"="Paraphyly"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.generalcategories.onlycompletelyaligned.allgenomes.png", device="png",dpi=300, width = 3, height = 3)

full.fractionation.status$M1.TimingNode<-full.fractionation.status$M1.TimingNode %>% 
  factor(levels = c("N0_basal","N0_Retained","N1_Tripsacum","N1_TdFL","N1_TdKS","N1_Zea","N2_Zn" ,"N2_ZdZm",
                    "N3_Zd","N3_Zm","N4_ZdGigi","N4_ZdMomo","N5_Zxvm",
                    "N5_Zh","N6_Zx","N6_Zvm","N8_Zv","N8_Zm","N7_ZxTIL18",
                    "N7_ZxTIL25","N9_ZvTIL11","N9_ZvTIL01",
                    "N10_Tropical","N10_Temperate", "N11_TropNoTzi8","N12_NCLines","N12_CMLandKLines","N13_NotCML333",
                    "N14_NotCML322","N15_CMLOnly","N15_CMLandKLines","N16_CML52CML69","N16_CML228andKLines","N17_NotKi11","N18_NotCML103",
                    "N19_NotM162W","N20_Ky21andOhLines","N20_FlintsBlinesMS71","N21_Flints","N21_BLinesandMS71","N22_NotB73","N23_NotKy21","N24_SweetNotHP301",
                    "ZmCML103","ZmCML228","ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmCML52","ZmCML69","ZmKi3","ZmNC350","ZmNC358","ZmTzi8",
                    "HP301",
                    "ZmB73","ZmB97","ZmKy21","ZmM162W","ZmMS71","ZmOh43","ZmOh7b","ZmP39","paraphyly"))
full.fractionation.status$M2.TimingNode<-full.fractionation.status$M2.TimingNode %>% 
  factor(levels = c("N0_basal","N0_Retained","N1_Tripsacum","N1_TdFL","N1_TdKS","N1_Zea","N2_Zn" ,"N2_ZdZm",
                    "N3_Zd","N3_Zm","N4_ZdGigi","N4_ZdMomo","N5_Zxvm",
                    "N5_Zh","N6_Zx","N6_Zvm","N8_Zv","N8_Zm","N7_ZxTIL18",
                    "N7_ZxTIL25","N9_ZvTIL11","N9_ZvTIL01","N10_Tropical","N10_Temperate", "N11_TropNoTzi8","N12_NCLines","N12_CMLandKLines","N13_NotCML333",
                    "N14_NotCML322","N15_CMLOnly","N15_CMLandKLines","N16_CML52CML69","N16_CML228andKLines","N17_NotKi11","N18_NotCML103",
                    "N19_NotM162W","N20_Ky21andOhLines","N20_FlintsBlinesMS71","N21_Flints","N21_BLinesandMS71","N22_NotB73","N23_NotKy21","N24_SweetNotHP301",
                    "ZmCML103","ZmCML228","ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmCML52","ZmCML69","ZmKi3","ZmNC350","ZmNC358","ZmTzi8",
                    "HP301",
                    "ZmB73","ZmB97","ZmKy21","ZmM162W","ZmMS71","ZmOh43","ZmOh7b","ZmP39","paraphyly"))

full.fractionation.status$M1.TimingNode.NoZnZd<-full.fractionation.status$M1.TimingNode.NoZnZd %>% 
  factor(levels = c("N0_basal","N0_Retained","N1_Tripsacum","N1_TdFL","N1_TdKS","N1_Zea","N3_Zd","N3_Zm","N4_ZdGigi","N4_ZdMomo","N5_Zxvm",
                    "N5_Zh","N6_Zx","N6_Zvm","N8_Zv","N8_Zm","N7_ZxTIL18",
                    "N7_ZxTIL25","N9_ZvTIL11","N9_ZvTIL01","N10_Tropical","N10_Temperate", "N11_TropNoTzi8","N12_NCLines","N12_CMLandKLines","N13_NotCML333",
                    "N14_NotCML322","N15_CMLOnly","N15_CMLandKLines","N16_CML52CML69","N16_CML228andKLines","N17_NotKi11","N18_NotCML103",
                    "N19_NotM162W","N20_Ky21andOhLines","N20_FlintsBlinesMS71","N21_Flints","N21_BLinesandMS71","N22_NotB73","N23_NotKy21","N24_SweetNotHP301",
                    "ZmCML103","ZmCML228","ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmCML52","ZmCML69","ZmKi3","ZmNC350","ZmNC358","ZmTzi8",
                    "HP301",
                    "ZmB73","ZmB97","ZmKy21","ZmM162W","ZmMS71","ZmOh43","ZmOh7b","ZmP39","paraphyly"))
full.fractionation.status$M2.TimingNode.NoZnZd<-full.fractionation.status$M2.TimingNode.NoZnZd %>% 
  factor(levels = c("N0_basal","N0_Retained","N1_Tripsacum","N1_TdFL","N1_TdKS","N1_Zea","N3_Zd","N3_Zm","N4_ZdGigi","N4_ZdMomo","N5_Zxvm",
                    "N5_Zh","N6_Zx","N6_Zvm","N8_Zv","N8_Zm","N7_ZxTIL18",
                    "N7_ZxTIL25","N9_ZvTIL11","N9_ZvTIL01","N10_Tropical","N10_Temperate", "N11_TropNoTzi8","N12_NCLines","N12_CMLandKLines","N13_NotCML333",
                    "N14_NotCML322","N15_CMLOnly","N15_CMLandKLines","N16_CML52CML69","N16_CML228andKLines","N17_NotKi11","N18_NotCML103",
                    "N19_NotM162W","N20_Ky21andOhLines","N20_FlintsBlinesMS71","N21_Flints","N21_BLinesandMS71","N22_NotB73","N23_NotKy21","N24_SweetNotHP301",
                    "ZmCML103","ZmCML228","ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmCML52","ZmCML69","ZmKi3","ZmNC350","ZmNC358","ZmTzi8",
                    "HP301",
                    "ZmB73","ZmB97","ZmKy21","ZmM162W","ZmMS71","ZmOh43","ZmOh7b","ZmP39","paraphyly"))

tibble(Node= c(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode),
       M = c(rep("M1",nrow(full.fractionation.status)),rep("M2",nrow(full.fractionation.status)))) %>%
  filter(!is.na(Node)) %>%
  ggplot(aes(y=Node))+
  geom_bar(aes(fill = M),position=position_dodge())+
  scale_fill_manual(values=subgenome_colors)+
  theme_minimal()+xlab("Count by Exon")+
  ggtitle("Timing using Exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.nodes.allgenomes.png", dpi=300, height=5, width=3, device="png")

tibble(Timing = c(full.fractionation.status$M1.Timing.NoZnZd, full.fractionation.status$M2.Timing.NoZnZd),
       M = c(rep("M1",nrow(full.fractionation.status)),rep("M2",nrow(full.fractionation.status)))) %>%
  filter(!is.na(Timing))%>%
  ggplot(aes(y=Timing))+
  geom_bar(aes(fill = M),position=position_dodge())+
  scale_fill_manual(values=subgenome_colors)+
  theme_minimal()+xlab("Count by Exon")+
  ggtitle("Timing using Exon without Zn")+
  scale_y_discrete(labels = c("basal" = "Basal", "genus"="Genus","genus_sister"="Genus, sister","species"="Species","subspecies_private"="Subspecies, private","accession_sister"="Accession, sister","private"="Private","CompletelyRetained"="Completely Retained","paraphyly"="Paraphyly","Some_Unaligned"="Some Unaligned","Completely_Unaligned"="Completely Unaligned"))

ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.generalcategories.noZnZd.png", device="png",dpi=300, width = 3, height = 3)
tibble(Timing = c(full.fractionation.status$M1.Timing.NoZnZd, full.fractionation.status$M2.Timing.NoZnZd),
       M = c(rep("M1",nrow(full.fractionation.status)),rep("M2",nrow(full.fractionation.status)))) %>%
  filter(!is.na(Timing) & !Timing %in% c("Completely_Unaligned","Some_Unaligned"))%>%
  ggplot(aes(y=Timing))+
  geom_bar(aes(fill = M),position=position_dodge())+
  scale_fill_manual(values=subgenome_colors)+
  theme_minimal()+xlab("Count by Exon")+
  ggtitle("Timing using Exon without Zn")+
  scale_y_discrete(labels = c("basal" = "Basal", "genus"="Genus","genus_sister"="Genus, sister","species"="Species","subspecies_private"="Subspecies, private","accession_sister"="Accession, sister","private"="Private","CompletelyRetained"="Completely Retained","paraphyly"="Paraphyly","Some_Unaligned"="Some Unaligned","Completely_Unaligned"="Completely Unaligned"))

ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.generalcategories.onlycompletelyaligned.noZnZd.png", device="png",dpi=300, width = 3, height = 3)

tibble(Node= c(full.fractionation.status$M1.TimingNode.NoZnZd, full.fractionation.status$M2.TimingNode.NoZnZd),
       M = c(rep("M1",nrow(full.fractionation.status)),rep("M2",nrow(full.fractionation.status)))) %>%
  filter(!is.na(Node)) %>%
  ggplot(aes(y=Node))+
  geom_bar(aes(fill = M),position=position_dodge())+
  scale_fill_manual(values=subgenome_colors)+
  theme_minimal()+xlab("Count by Exon")+
  ggtitle("Timing using Exon without Zn")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.nodes.noZnZd.png", dpi=300, height=5, width=3, device="png")

#(double check that all the way up ran after adding in the complete retention)
#how many exons are in N0 vs any other?
filter(full.fractionation.status, M1.TimingNode == "N0_basal") %>% nrow() #12553
filter(full.fractionation.status, M1.TimingNode != "N0_basal") %>% nrow() #32251
filter(full.fractionation.status, M2.TimingNode == "N0_basal") %>% nrow() #11938
filter(full.fractionation.status, M2.TimingNode != "N0_basal") %>% nrow() #14880

filter(full.fractionation.status, M1.TimingNode.NoZnZd == "N0_basal") %>% nrow() #16620
filter(full.fractionation.status, M1.TimingNode.NoZnZd != "N0_basal") %>% nrow() #42544
filter(full.fractionation.status, M2.TimingNode.NoZnZd == "N0_basal") %>% nrow() #20033
filter(full.fractionation.status, M2.TimingNode.NoZnZd != "N0_basal") %>% nrow() #24288

#How does removing Zn change the timing estimates? Change the amount of paraphyly?

filter(full.fractionation.status, M1.Timing == "paraphyly") %>% nrow() #4826
filter(full.fractionation.status, M1.Timing != "paraphyly") %>% nrow() #64440
filter(full.fractionation.status, M2.Timing == "paraphyly") %>% nrow() #3326
filter(full.fractionation.status, M2.Timing != "paraphyly") %>% nrow() #65933

filter(full.fractionation.status, M1.Timing.NoZnZd == "paraphyly") %>% nrow() #6231
filter(full.fractionation.status, M1.Timing.NoZnZd != "paraphyly") %>% nrow() #63030
filter(full.fractionation.status, M2.Timing.NoZnZd == "paraphyly") %>% nrow() #5245
filter(full.fractionation.status, M2.Timing.NoZnZd != "paraphyly") %>% nrow() #64011

#IDK if this will get me the answer I'm looking for...
#But suggests only a handful of exons change their timing estimates when excluding Zn
filter(full.fractionation.status, as.character(M1.TimingNode) != as.character(M1.TimingNode.NoZnZd)) %>% nrow() #852
filter(full.fractionation.status, as.character(M2.TimingNode) != as.character(M2.TimingNode.NoZnZd)) %>% nrow() #285

#Is there a difference in the rate of fractionation between M1 and M2?
table(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode) %>% view()

fisher.test(table(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode), simulate.p.value = T) 
#Fisher's Exact Test for Count Data with simulated p-value (based on 2000 replicates)

#data:  table(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode)
#p-value = 0.0004998
#alternative hypothesis: two.sided

fisher.test(table(full.fractionation.status$M1.TimingNode.NoZnZd, full.fractionation.status$M2.TimingNode.NoZnZd), simulate.p.value = T) 
#data:  table(full.fractionation.status$M1.TimingNode.NoZnZd, full.fractionation.status$M2.TimingNode.NoZnZd)
#p-value = 0.0004998
#alternative hypothesis: two.sided

#Why is there a N11_Tropical???

table(full.fractionation.status$M1.TimingNode.NoZnZd, full.fractionation.status$M2.TimingNode.NoZnZd) %>% as.data.frame()%>%
  ggplot(aes(x=Var1, y=Var2, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(option="A", values = scales::rescale(x = c(0,1:6 %o% 10^(1:3)), to = c(0,1)), direction = -1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  xlab("M1")+ylab("M2")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.nodes.sharing.noZnZd.png", device="png",dpi=300,height = 8, width = 8.5)

table(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode) %>% as.data.frame()%>%
  ggplot(aes(x=Var1, y=Var2, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(option="A", values = scales::rescale(x = c(0,1:6 %o% 10^(1:3)), to = c(0,1)), direction = -1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  xlab("M1")+ylab("M2")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.nodes.sharing.allGenomes.png", device="png",dpi=300,height = 8, width = 8.5)

#Changing to percentage? 
table(full.fractionation.status$M1.TimingNode.NoZnZd, full.fractionation.status$M2.TimingNode.NoZnZd) %>% as.data.frame() %>% pull(Freq) %>% sum()
#denominator is 38554
table(full.fractionation.status$M1.TimingNode.NoZnZd, full.fractionation.status$M2.TimingNode.NoZnZd) %>% 
  as.data.frame() %>%
  mutate(Percent = (Freq / 38554)*100) %>% 
  ggplot(aes(x=Var1, y=Var2, fill=Percent))+
  geom_tile()+
  scale_fill_viridis_c(option="A", values = scales::rescale(x = c(0,1:6 %o% 10^(1:3)), to = c(0,1)), direction = -1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  xlab("M1")+ylab("M2")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.nodes.sharing.noZnZd.percent.png", device="png",dpi=300,height = 8, width = 8.5)


table(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode) %>% as.data.frame() %>% pull(Freq) %>% sum()
#denominator = 16603
table(full.fractionation.status$M1.TimingNode, full.fractionation.status$M2.TimingNode) %>% as.data.frame()%>%
  mutate(Percent = (Freq / 16603)*100) %>% 
  ggplot(aes(x=Var1, y=Var2, fill=Percent))+
  geom_tile()+
  scale_fill_viridis_c(option="A", values = scales::rescale(x = c(0,1:6 %o% 10^(1:3)), to = c(0,1)), direction = -1)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90))+
  xlab("M1")+ylab("M2")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/timing.nodes.sharing.allGenomes.percent..png", device="png",dpi=300,height = 8, width = 8.5)

#How is paraphyly split across lineages?

paraphyly.M1<-filter(full.fractionation.status, M1.Timing == "paraphyly")
#get rid of extra columns
paraphyly.M1<-paraphyly.M1 %>% select(-contains("M2"))

paraphyly.M1 %>% select(c(ends_with(".M1"))) %>% .[1,] %>% as.data.frame() %>% select(where(~ .x == 1)) %>% colnames()%>% str_remove(".M1")%>% paste(collapse =":")

paraphyly.M1$paraphyly_type<-NA
for(i in 1:nrow(paraphyly.M1)){
  paraphyly.M1$paraphyly_type[i]<-paraphyly.M1[i,] %>% as.data.frame()%>% select(c(ends_with(".M1"))) %>% select(where(~ .x == 1))%>% colnames()%>% str_remove(".M1")%>% paste(collapse =":")
}

paraphyly.M2<-filter(full.fractionation.status, M2.Timing == "paraphyly")
paraphyly.M2<-paraphyly.M2 %>% select(-contains("M1"))
paraphyly.M2$paraphyly_type<-NA
for(i in 1:nrow(paraphyly.M2)){
  paraphyly.M2$paraphyly_type[i]<-paraphyly.M2[i,] %>% as.data.frame()%>% select(c(ends_with(".M2"))) %>% select(where(~ .x == 1))%>% colnames()%>% str_remove(".M2")%>% paste(collapse =":")
}

paraphyly.M1.NoZnZd<-filter(full.fractionation.status, M1.Timing.NoZnZd == "paraphyly")
paraphyly.M1.NoZnZd<-paraphyly.M1.NoZnZd %>% select(-c(contains("M2"),contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))
paraphyly.M1.NoZnZd$paraphyly_type<-NA
for(i in 1:nrow(paraphyly.M1.NoZnZd)){
  paraphyly.M1.NoZnZd$paraphyly_type[i]<-paraphyly.M1.NoZnZd[i,] %>% as.data.frame() %>% select(c(ends_with(".M1"))) %>% select(where(~.x == 1)) %>% colnames() %>% str_remove(".M1") %>% paste(collapse = ":")
}
paraphyly.M2.NoZnZd<-filter(full.fractionation.status, M2.Timing.NoZnZd == "paraphyly")
paraphyly.M2.NoZnZd$paraphyly_type<-NA
paraphyly.M2.NoZnZd<-paraphyly.M2.NoZnZd %>% select(-c(contains("M1"),contains("Zn"),starts_with("ZdGigi."),starts_with("ZdMomo.")))
for(i in 1:nrow(paraphyly.M2.NoZnZd)){
  paraphyly.M2.NoZnZd$paraphyly_type[i]<-paraphyly.M2.NoZnZd[i,] %>% as.data.frame() %>% select(c(ends_with(".M2"))) %>% select(where(~.x == 1)) %>% colnames() %>% str_remove(".M2") %>% paste(collapse = ":")
}

#Number of types of paraphyly:
paraphyly.M1$paraphyly_type %>% unique() %>% length() #1847 paraphyly types, M1 all genomes
paraphyly.M1.NoZnZd$paraphyly_type %>% unique() %>% length() #2354 paraphyly types, M1 no ZnZd (how'd it get worse? NA's in Zn and Zd probably)
paraphyly.M2$paraphyly_type %>% unique() %>% length() #1267 paraphyly types, M2 all genomes
paraphyly.M2.NoZnZd$paraphyly_type %>% unique() %>% length() #1863 paraphyly types, M2 noZnZd (a little bit better, but not much!)


paraphyly.M1.summary<-paraphyly.M1 %>% group_by(paraphyly_type) %>% count()
paraphyly.M1.NoZnZd.summary<-paraphyly.M1.NoZnZd %>% group_by(paraphyly_type) %>% count()
paraphyly.M2.summary<-paraphyly.M2 %>% group_by(paraphyly_type) %>% count()
paraphyly.M2.NoZnZd.summary<-paraphyly.M2.NoZnZd %>% group_by(paraphyly_type) %>% count()

summary(paraphyly.M1.summary$n)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   2.613   1.000 165.000 
summary(paraphyly.M2.summary$n)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   2.625   2.000 139.000 

summary(paraphyly.M1.NoZnZd.summary$n)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   2.647   1.000 204.000 
summary(paraphyly.M2.NoZnZd.summary$n)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   2.815   1.000 197.000 

#Removing Zn and Zd doesn't drastically change the average number of paraphyly patterns
#(M1 max goes from 126 --> 135, M2 max goes from 59--> 83), likely due to fewer NAs tossing exons out

#What's the number of genomes involved in a given paraphyly pattern?
paraphyly.M1.summary<-paraphyly.M1.summary %>% mutate(NumGenomesInvolved = str_split(paraphyly_type, ":",simplify = T)%>% length())
paraphyly.M2.summary<-paraphyly.M2.summary %>% mutate(NumGenomesInvolved = str_split(paraphyly_type, ":",simplify = T)%>% length())
paraphyly.M1.NoZnZd.summary<-paraphyly.M1.NoZnZd.summary %>% mutate(NumGenomesInvolved = str_split(paraphyly_type, ":",simplify = T)%>% length())
paraphyly.M2.NoZnZd.summary<-paraphyly.M2.NoZnZd.summary %>% mutate(NumGenomesInvolved = str_split(paraphyly_type, ":",simplify = T)%>% length())

summary(paraphyly.M1.summary$NumGenomesInvolved)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.00    6.00   16.00   18.09   30.00   38.00  
summary(paraphyly.M2.summary$NumGenomesInvolved)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    8.00   20.00   19.64   32.00   38.00 

summary(paraphyly.M1.NoZnZd.summary$NumGenomesInvolved)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#2.00    6.00   16.00   16.72   28.00   34.00 
summary(paraphyly.M2.NoZnZd.summary$NumGenomesInvolved)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    7.00   18.00   17.33   28.00   34.00 

#Removing ZnZd changes the max (because 3 less genomes possible [38 --> 34])
#The mean and median also drop by about 1, M2 becomes more like M1

detect_genome_in_paraphyly_type<-function(df){
  new_df<-mutate(df, TdFL = str_detect(paraphyly_type, "TdFL"),
                 TdKS = str_detect(paraphyly_type, "TdKS"),
                 ZdGigi_4to1 = str_detect(paraphyly_type, "ZdGigi_4to1"),
                 ZdMomo_4to1 = str_detect(paraphyly_type, "ZdMomo_4to1"),
                 ZdGigi = str_detect(paraphyly_type, "ZdGigi\\b"),
                 ZdMomo = str_detect(paraphyly_type, "ZdMomo\\b"),
                 ZhRIMHU001 = str_detect(paraphyly_type, "ZhRIMHU001"),
                 ZmB73 = str_detect(paraphyly_type, "ZmB73"),
                 ZmB97 = str_detect(paraphyly_type, "ZmB97"),
                 ZmCML103 = str_detect(paraphyly_type, "ZmCML103"),
                 ZmCML228 = str_detect(paraphyly_type, "ZmCML228"),
                 ZmCML247 = str_detect(paraphyly_type, "ZmCML247"),
                 ZmCML277 = str_detect(paraphyly_type, "ZmCML277"),
                 ZmCML322 = str_detect(paraphyly_type, "ZmCML322"),
                 ZmCML333 = str_detect(paraphyly_type, "ZmCML333"),
                 ZmCML52 = str_detect(paraphyly_type, "ZmCML52"),
                 ZmCML69 = str_detect(paraphyly_type, "ZmCML69"),
                 ZmHP301 = str_detect(paraphyly_type, "ZmHP301"),
                 ZmIL14H = str_detect(paraphyly_type, "ZmIL14H"),
                 ZmKi11 = str_detect(paraphyly_type, "ZmKi11"),
                 ZmKi3 = str_detect(paraphyly_type, "ZmKi3"),
                 ZmKy21 = str_detect(paraphyly_type, "ZmKy21"),
                 ZmM162W = str_detect(paraphyly_type, "ZmM162W"),
                 
                
                 ZmMS71 = str_detect(paraphyly_type, "ZmMS71"),
                 ZmNC350 = str_detect(paraphyly_type, "ZmNC350"),
                 ZmNC358 = str_detect(paraphyly_type, "ZmNC358"),
                 ZmOh43 = str_detect(paraphyly_type, "ZmOh43"),
                 ZmOh7b = str_detect(paraphyly_type, "ZmOh7b"),
                 ZmP39 = str_detect(paraphyly_type, "ZmP39"),
                 
                 ZmTzi8 = str_detect(paraphyly_type, "ZmTzi8"),
                 ZnPI615697 = str_detect(paraphyly_type, "ZnPI615697\\b"),
                 ZnPI615697_4to1 = str_detect(paraphyly_type, "ZnPI615697_4to1"),
                 ZvTIL01 = str_detect(paraphyly_type, "ZvTIL01"),
                 ZvTIL11 = str_detect(paraphyly_type, "ZvTIL11"),
                 ZxTIL18 = str_detect(paraphyly_type, "ZxTIL18"),
                 ZxTIL25 = str_detect(paraphyly_type, "ZxTIL25"))
  return(new_df)}
paraphyly.M1.summary<-detect_genome_in_paraphyly_type(paraphyly.M1.summary)
paraphyly.M2.summary<-detect_genome_in_paraphyly_type(paraphyly.M2.summary)

paraphyly.M1.NoZnZd.summary<-detect_genome_in_paraphyly_type(paraphyly.M1.NoZnZd.summary)%>% select(-c("ZdGigi","ZdMomo",contains("Zn")))
paraphyly.M2.NoZnZd.summary<-detect_genome_in_paraphyly_type(paraphyly.M2.NoZnZd.summary)%>% select(-c("ZdGigi","ZdMomo",contains("Zn")))

lapply(paraphyly.M1.summary[,4:ncol(paraphyly.M1.summary)], sum) %>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M1")+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.AllGenomes.png", dpi=300, device="png",width=5, height=4)

lapply(paraphyly.M2.summary[,4:ncol(paraphyly.M2.summary)], sum) %>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M2")+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.AllGenomes.png", dpi=300, device="png",width=5,height=4)

lapply(paraphyly.M1.NoZnZd.summary[,4:ncol(paraphyly.M1.NoZnZd.summary)], sum) %>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M1")+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.NoZnZd.png", dpi=300, device="png",width=5, height=4)

lapply(paraphyly.M2.NoZnZd.summary[,4:ncol(paraphyly.M2.NoZnZd.summary)], sum) %>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M2")+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.NoZnZd.png", dpi=300, device="png",width=5,height=4)

full.paraphyly<-full_join(x=paraphyly.M1.summary, y=paraphyly.M2.summary, by="paraphyly_type", suffix = c(".M1",".M2"))
full.paraphyly.NoZnZd<-full_join(x=paraphyly.M1.NoZnZd.summary, y=paraphyly.M2.NoZnZd.summary, by="paraphyly_type", suffix = c(".M1",".M2"))

#are there any patterns that are present in M1 or M2 but not the other?
nrow(full.paraphyly) #3049
filter(full.paraphyly, is.na(n.M1) | is.na(n.M2)) %>% nrow() #filters from 3049 patterns to 2984 (so x patterns shared between M1 and M2)
filter(full.paraphyly,  is.na(n.M2)) %>% nrow() #is there in M1 but not M2 #1782
filter(full.paraphyly, is.na(n.M1)) %>% nrow() #is there in M2 but not M1 #1202

nrow(full.paraphyly.NoZnZd) #4127
filter(full.paraphyly.NoZnZd, is.na(n.M1) | is.na(n.M2)) %>% nrow() #filters from 4127 patterns to 4037 (so x patterns shared between M1 and M2)
filter(full.paraphyly.NoZnZd,  is.na(n.M2))  %>% nrow()#is there in M1 but not M2 #2264
filter(full.paraphyly.NoZnZd, is.na(n.M1))   %>% nrow()#is there in M2 but not M1 #1773

#What are the patterns that are above the 3rd quartile in terms of frequency? Still may be too large to look at...
paraphyly.M1.NoZnZd.summary$NumGenomesInvolved %>% summary()
filter(paraphyly.M1.NoZnZd.summary, n >=25)%>% select(NumGenomesInvolved) %>% summary() 
paraphyly.M2.NoZnZd.summary$NumGenomesInvolved %>% summary()
filter(paraphyly.M2.NoZnZd.summary, n >=25)%>% select(NumGenomesInvolved) %>% summary() 

#M1 all genomes == 3
filter(paraphyly.M1.summary, n >= 3) 
filter(paraphyly.M1.summary, n >= 3) %>% select(NumGenomesInvolved) %>% summary() #average number of genomes involved goes up!
#M2 all genomes == 4
filter(paraphyly.M2.summary, n >= 4) 
filter(paraphyly.M2.summary, n >= 4) %>% select(NumGenomesInvolved) %>% summary() #average number of genomes involved goes up!
#M1 No ZnZd == 3
filter(paraphyly.M1.NoZnZd.summary, n >= 3) 
filter(paraphyly.M1.NoZnZd.summary, n >= 3) %>% select(NumGenomesInvolved) %>% summary() #average number of genomes involved goes up!
#M2 No ZnZd == 4
filter(paraphyly.M2.NoZnZd.summary, n >= 4) 
filter(paraphyly.M2.NoZnZd.summary, n >= 4) %>% select(NumGenomesInvolved) %>% summary() #average number of genomes involved goes up!

#Let's see which patterns occur more than 100 times?
#M1 all genomes == 3
filter(paraphyly.M1.summary, n >= 100) %>% select(paraphyly_type, n, NumGenomesInvolved) %>% arrange(-n)
filter(paraphyly.M1.summary, n >= 100) %>%  arrange(-n) %>% view()
filter(paraphyly.M1.summary, n >= 100) %>% ungroup() %>% select(-c(paraphyly_type, n, NumGenomesInvolved)) %>% lapply(sum)
#M2 all genomes == 4
filter(paraphyly.M2.summary, n >= 100) %>% select(paraphyly_type, n, NumGenomesInvolved) %>% arrange(-n)
filter(paraphyly.M2.summary, n >= 100) %>% ungroup() %>% select(-c(paraphyly_type, n, NumGenomesInvolved)) %>% lapply(sum)
#M1 No ZnZd == 3
filter(paraphyly.M1.NoZnZd.summary, n >= 100) %>% select(paraphyly_type, n, NumGenomesInvolved) %>% arrange(-n)
filter(paraphyly.M1.NoZnZd.summary, n >= 100) %>% ungroup() %>% select(-c(paraphyly_type, n, NumGenomesInvolved)) %>% lapply(sum)

#M2 No ZnZd == 4
filter(paraphyly.M2.NoZnZd.summary, n >= 100) %>% select(paraphyly_type, n, NumGenomesInvolved) %>% arrange(-n)
filter(paraphyly.M2.NoZnZd.summary, n >= 100) %>% ungroup() %>% select(-c(paraphyly_type, n, NumGenomesInvolved)) %>% lapply(sum)

filter(paraphyly.M1.NoZnZd.summary, n >= 100) %>% arrange(-n)%>%
  ggplot(aes(y=NumGenomesInvolved))+
  geom_point(aes(x=n))#+
#theme(axis.text.y = element_blank())+ylab("")

ggplot(full.paraphyly)+
  geom_point(aes(x=NumGenomesInvolved.M1, y=n.M1), color=subgenome_colors[1], alpha=0.5)+
  geom_point(aes(x=NumGenomesInvolved.M2, y=n.M2), color=subgenome_colors[2], alpha=0.5)+
  scale_color_manual(values = subgenome_colors)+
  labs(color="M")

ggplot(full.paraphyly)+
  geom_point(aes(x=NumGenomesInvolved.M1, y=n.M1), color=subgenome_colors[1], alpha=0.5)+
  labs(color="M",x="Number of Genomes Involved", y="Count of Genes per Pattern")+
  theme_minimal()+
  ggtitle("M1 Paraphyly Patterns, All Genomes")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.PatternCnts.AllGenomes.png", dpi=300, width=3, height=3)
ggplot(full.paraphyly)+
  geom_point(aes(x=NumGenomesInvolved.M2, y=n.M2), color=subgenome_colors[2], alpha=0.5)+
  labs(color="M",x="Number of Genomes Involved", y="Count of Genes per Pattern")+
  theme_minimal()+
  ggtitle("M2 Paraphyly Patterns, All Genomes")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.PatternCnts.AllGenomes.png", dpi=300, width=3, height=3)

ggplot(full.paraphyly.NoZnZd)+
  geom_point(aes(x=NumGenomesInvolved.M1, y=n.M1), color=subgenome_colors[1], alpha=0.5)+
  labs(color="M",x="Number of Genomes Involved", y="Count of Genes per Pattern")+
  theme_minimal()+
  ggtitle("M1 Paraphyly Patterns, No Zn")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.PatternCnts.NoZnZd.png", dpi=300, width=3, height=3)
ggplot(full.paraphyly.NoZnZd)+
  geom_point(aes(x=NumGenomesInvolved.M2, y=n.M2), color=subgenome_colors[2], alpha=0.5)+
  labs(color="M",x="Number of Genomes Involved", y="Count of Genes per Pattern")+
  theme_minimal()+
  ggtitle("M2 Paraphyly Patterns, No Zn")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.PatternCnts.NoZnZd.png", dpi=300, width=3, height=3)

filter(paraphyly.M2.summary, NumGenomesInvolved < 5) %>% .[,4:ncol(paraphyly.M2.summary)] %>% lapply(sum)%>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M2")+ylab("")
ggsave('/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.NumGenomeInvolvedlt5.AllGenomes.png',dpi=300, width=5,height=4)

filter(paraphyly.M1.summary, NumGenomesInvolved < 5) %>% .[,4:ncol(paraphyly.M1.summary)] %>% lapply(sum)%>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M1")+ylab("")
ggsave('/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.NumGenomeInvolvedlt5.AllGenomes.png',dpi=300, width=5,height=4)

filter(paraphyly.M2.NoZnZd.summary, NumGenomesInvolved < 5) %>% .[,4:ncol(paraphyly.M2.NoZnZd.summary)] %>% lapply(sum)%>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M2")+ylab("")
ggsave('/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.NumGenomeInvolvedlt5.noZnZd.png',dpi=300, width=5,height=4)

filter(paraphyly.M1.NoZnZd.summary, NumGenomesInvolved < 5) %>% .[,4:ncol(paraphyly.M1.NoZnZd.summary)] %>% lapply(sum)%>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  xlab("Number of Paraphyly Patterns Including Each Genome, M1")+ylab("")
ggsave('/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.NumGenomeInvolvedlt5.NoZnZd.png',dpi=300, width=5,height=4)

#What about taking the patterns that involve 25 genes? Exons? or more
filter(paraphyly.M1.NoZnZd.summary, n >= 25) %>% .[,4:ncol(paraphyly.M1.NoZnZd.summary)] %>% lapply(sum)%>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("Number of Paraphyly Patterns Including Each Genome, M1")+ylab("")
ggsave('/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M1.GenomesInvolved.Ngt25.NoZnZd.png',dpi=300, width=5,height=4)

filter(paraphyly.M2.NoZnZd.summary, n >= 25) %>% .[,4:ncol(paraphyly.M2.NoZnZd.summary)] %>% lapply(sum)%>% as_tibble() %>%
  pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")))%>%
  ggplot(aes(x=NumParaphylyPatterns,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("Number of Paraphyly Patterns Including Each Genome, M2")+ylab("")
ggsave('/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/paraphyly.M2.GenomesInvolved.Ngt25.NoZnZd.png',dpi=300, width=5,height=4)

#what if the cut off is 30 and plot the inverse?
filter(paraphyly.M2.NoZnZd.summary, n >= 30) %>% nrow() #14
filter(paraphyly.M2.NoZnZd.summary, n >= 30) %>% .[,4:ncol(paraphyly.M2.NoZnZd.summary)] %>% lapply(sum)%>% as_tibble() %>%pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")),
         Inverse_involvement = 14-NumParaphylyPatterns)%>%
  ggplot(aes(x=Inverse_involvement,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("Count of Paraphyly Patterns Excluding Each Genome, M2")+ylab("")+ggtitle("Paraphyly Pattern Involves 30+ Genomes")
  
filter(paraphyly.M1.NoZnZd.summary, n >= 30) %>% nrow() #16
filter(paraphyly.M1.NoZnZd.summary, n >= 30) %>% .[,4:ncol(paraphyly.M2.NoZnZd.summary)] %>% lapply(sum)%>% as_tibble() %>%pivot_longer(cols = everything(),names_to = "Genome",values_to = "NumParaphylyPatterns")%>%
  mutate(Genome = Genome %>% factor(levels=c("ZmB73" ,"ZmB97", "ZmIL14H","ZmKy21" ,"ZmM162W" ,"ZmMS71" ,"ZmOh43","ZmOh7b","ZmP39" ,
                                             "ZmHP301" ,
                                             "ZmCML103" ,"ZmCML228" ,"ZmCML247" ,"ZmCML277" ,"ZmCML322" ,"ZmCML333" ,"ZmCML52" ,"ZmCML69" ,"ZmKi11" ,"ZmKi3" ,"ZmNC350" ,"ZmNC358" ,"ZmTzi8" ,
                                             "ZvTIL11" ,"ZvTIL01" ,"ZxTIL18" ,"ZxTIL25","ZhRIMHU001" ,"ZdGigi","ZdMomo","ZnPI615697","ZdGigi_4to1","ZdMomo_4to1","ZnPI615697_4to1","TdFL","TdKS")),
         Inverse_involvement = 16-NumParaphylyPatterns)%>%
  ggplot(aes(x=Inverse_involvement,y=Genome))+
  geom_bar(stat="identity",aes(fill=Genome))+
  scale_fill_manual(values=genome_colors)+
  theme_minimal()+
  theme(legend.position = "none")+
  xlab("Count of Paraphyly Patterns Excluding Each Genome, M1")+ylab("")+ggtitle("Paraphyly Pattern Involves 30+ Genomes")

####ALSO WHAT IS THE NUMBER OF M1/M2 EXONS THAT GO AT EACH NODE?####
Timing.M1<-full.fractionation.status %>% select(c(Gene_ID, ID, contains("M1.Timing")))
Timing.M2<-full.fractionation.status %>% select(c(Gene_ID, ID, contains("M2.Timing")))

Timing.M1 %>% group_by(M1.TimingNode) %>% filter(!is.na(M1.TimingNode)) %>% count() %>% select(n) %>% pull() %>% sum() #44804
Timing.M2 %>% group_by(M2.TimingNode) %>% filter(!is.na(M2.TimingNode)) %>% count() %>% select(n) %>% pull() %>% sum() #26818

Timing.M1 %>% group_by(M1.TimingNode) %>% count() %>% mutate(PercentofAllExons = (n/69269)*100, 
                                                             PercentofCompletelyAlignedM1Exons = (n/44804)*100)
Timing.M2 %>% group_by(M2.TimingNode) %>% count() %>% mutate(PercentofAllExons = (n/69269)*100, 
                                                             PercentofCompletelyAlignedM2Exons = (n/26818)*100)
Timing.M1 %>% filter(M1.Timing %in% c("Some_Unaligned","Completely_Unaligned")) %>% group_by(M1.Timing) %>% count() %>% mutate(percent = (n/69269)*100)
Timing.M2 %>% filter(M2.Timing %in% c("Some_Unaligned","Completely_Unaligned")) %>% group_by(M2.Timing) %>% count()%>% mutate(percent = (n/69269)*100)


Timing.M1 %>% group_by(M1.TimingNode.NoZnZd) %>% filter(!is.na(M1.TimingNode.NoZnZd)) %>% count() %>% select(n) %>% pull() %>% sum() #59164
Timing.M2 %>% group_by(M2.TimingNode.NoZnZd) %>% filter(!is.na(M2.TimingNode.NoZnZd)) %>% count() %>% select(n) %>% pull() %>% sum() #44321

Timing.M1 %>% group_by(M1.TimingNode.NoZnZd) %>% count() %>% mutate(PercentofAllExons = (n/69269)*100, 
                                                                    PercentofCompletelyAlignedM1Exons = (n/59164)*100) %>%
  write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/Timing.M1.NoZnZd.tsv")
Timing.M2 %>% group_by(M2.TimingNode.NoZnZd) %>% count() %>% mutate(PercentofAllExons = (n/69269)*100, 
                                                                    PercentofCompletelyAlignedM2Exons = (n/44321)*100)%>%
  write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/Timing.M2.NoZnZd.tsv")
Timing.M1 %>% group_by(M1.TimingNode.NoZnZd) %>% count() %>% filter(!M1.TimingNode.NoZnZd %in% c("N0_basal","N0_Retained","N1_Tripsacum","N1_Zea","paraphyly") & !is.na(M1.TimingNode.NoZnZd)) %>% select(n) %>% pull() %>% sum()
Timing.M2 %>% group_by(M2.TimingNode.NoZnZd) %>% count() %>% filter(!M2.TimingNode.NoZnZd %in% c("N0_basal","N0_Retained","N1_Tripsacum","N1_Zea","paraphyly") & !is.na(M2.TimingNode.NoZnZd)) %>% select(n) %>% pull() %>% sum()

####VCF DELETION ANALYSIS####
deletion_stats<-tibble(RefChr = as.character(c(2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9)),
                       Subgenome =c("M1","M2","M1","M2","M1","M2","M1","M2","M1","M2","M1","M2","M1","M2","M1","M2"),
                       N_exonic_dels = c(8905,11388,11233,13341,8154,10785,2984,3414,6591,7487,4888,4681,3242,3846,6351,6468),#number of deletions within exons
                       N_exons_withdel = c(3663,5479,4088,6596,4012,5402,1128,1579,2704,3520,1980,2331,1270,1971,2362,3167),#number of exons with a deletion
                       N_single_dels = c(616,460,599,631,646,483,214,167,447,327,199,257,219,229,318,362), #number of deletions that are the only deletion within the exon
                       N_overlapping_dels = c(3047,5019,3489,5965,3366,3115,914,1412,2257,3193,1781,2074,1051,1742,2044,2805), #number of deletions that overlap each other
                       N_exactBoundaries =c(5978,7352,8041,8733,5626,7226,1849,2014,4521,4831,3277,2781,2021,2374,4471,4452), #number of deletions completely within the boundaries of an exon
                       N_inframe =c(1272,1823,1716,2184,1147,1691,453,615,1051,1100,670,674,415,610,942,989), #number of deletions that are in frame
                       N_frameshift =c(4706,5529,6325,6549,4479,5535,1396,1399,3470,3731,2607,2107,1606,1764,3529,3463),#number of frame shifting deletions
                       ) 
deletion_stats<-mutate(deletion_stats, Percent_single_dels = (N_single_dels/N_exonic_dels)*100,
                       Percent_overlapping_dels = (N_overlapping_dels/N_exonic_dels)*100,
                       Percent_exactBoundaries = (N_exactBoundaries/N_exonic_dels)*100,
                       Percent_inframe = (N_inframe/N_exactBoundaries)*100,
                       Percent_frameshift = (N_frameshift/N_exactBoundaries)*100)
ggplot(deletion_stats, aes(x=RefChr, y=Percent_single_dels))+
  geom_bar(stat = "identity")+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  ylab("Percent Deletions that are un-nested")

ggplot(deletion_stats, aes(x=RefChr, y=Percent_overlapping_dels))+
  geom_bar(stat = "identity")+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  ylab("Percent Deletions that are nested")

#stacked barchart
tibble(RefChr = rep(deletion_stats$RefChr,2),
       Subgenome = rep(deletion_stats$Subgenome,2),
       Type = c(rep("Single deletion", nrow(deletion_stats)),
                rep("Nested deletion", nrow(deletion_stats))),
       Percent = c(deletion_stats$Percent_single_dels, deletion_stats$Percent_overlapping_dels)) %>%
  ggplot(aes(x=RefChr, y=Percent, fill = Type))+
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  ylab("Percent")+ggtitle("Percent of deletions that are nested vs. single within an exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/vcf-dels.nestedvssingleton.png",device = "png",dpi=300, width = 4, height = 3)

#violin
tibble(RefChr = rep(deletion_stats$RefChr,2),
       Subgenome = rep(deletion_stats$Subgenome,2),
       Type = c(rep("Single deletion", nrow(deletion_stats)),
                rep("Nested deletion", nrow(deletion_stats))),
       Percent = c(deletion_stats$Percent_single_dels, deletion_stats$Percent_overlapping_dels)) %>%
  ggplot(aes(x=Subgenome, y=Percent, fill = Subgenome))+
  geom_violin()+
  geom_jitter(height=0,width = 0.1)+
  facet_wrap(vars(Type))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_fill_manual(values=subgenome_colors)+
  ylab("Percent")+ggtitle("Types of Deletions")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/vcf-dels.typesofdels.violin.png",device = "png",dpi=300, width = 4, height = 3)


ggplot(deletion_stats, aes(x=Subgenome, y=Percent_exactBoundaries,fill=Subgenome))+
  geom_violin()+
  geom_jitter(height = 0, width = 0.1)+
  scale_fill_manual(values = subgenome_colors)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Deletions contained within an exon")+ylab("Percent")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/vcf-dels.exactdeletions.violin.png",device = "png",dpi=300, width = 4, height = 4)

ggplot(deletion_stats, aes(x=RefChr, y=Percent_inframe))+
  geom_bar(stat = "identity")+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  ylab("Percent Deletions that are inframe")

ggplot(deletion_stats, aes(x=RefChr, y=Percent_frameshift))+
  geom_bar(stat = "identity")+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  ylab("Percent Deletions that are frame shifting")

#stacked barchart
tibble(RefChr = rep(deletion_stats$RefChr,2),
       Subgenome = rep(deletion_stats$Subgenome,2),
       Type = c(rep("In Frame", nrow(deletion_stats)),
                rep("Frame Shift", nrow(deletion_stats))),
       Percent = c(deletion_stats$Percent_inframe, deletion_stats$Percent_frameshift)) %>%
  ggplot(aes(x=RefChr, y=Percent, fill = Type))+
  geom_bar(stat = "identity", position = "stack")+
  facet_wrap(vars(Subgenome))+
  theme_bw()+
  ylab("Percent")+ggtitle("Frame shifting deletions within exons")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/vcf-dels.frameshift.png",device = "png",dpi=300, width = 4, height = 3)

#violin plot
tibble(RefChr = rep(deletion_stats$RefChr,1),
       Subgenome = rep(deletion_stats$Subgenome,1),
       Type = c(rep("Frame Shift", nrow(deletion_stats))),
       Percent = c(deletion_stats$Percent_frameshift)) %>%
  ggplot(aes(x=Subgenome, y=Percent, fill = Subgenome))+
  geom_violin()+
  geom_jitter(height = 0, width = 0.1)+
  scale_fill_manual(values = subgenome_colors)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Frame shifting deletions within an exon")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/vcf-dels.frameshiftOnly.violin.png",device = "png",dpi=300, width = 4, height = 4)

#put them all together
tibble(RefChr = rep(deletion_stats$RefChr,4),
       Subgenome = rep(deletion_stats$Subgenome,4),
       Type = c(rep("Single deletion", nrow(deletion_stats)),
                rep("Nested deletion", nrow(deletion_stats)),
                rep("Contained w/in Single Exon",nrow(deletion_stats)),
                rep("Frame Shift", nrow(deletion_stats))),
       Percent = c(deletion_stats$Percent_single_dels, 
                   deletion_stats$Percent_overlapping_dels,
                   deletion_stats$Percent_exactBoundaries,
                   deletion_stats$Percent_frameshift)) %>%
  ggplot(aes(x=Subgenome, y=Percent, fill = Subgenome))+
  geom_violin()+
  geom_jitter(height = 0, width = 0.1)+
  facet_wrap(vars(Type))+
  scale_fill_manual(values = subgenome_colors)+
  theme_bw()+
  theme(legend.position = "none")+
  ggtitle("Types of Deletions")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/vcf-dels.TypeOfDeletion.violin.png",device = "png",dpi=300, width = 4, height = 4)



