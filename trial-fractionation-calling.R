#To look at deletion calls as fractionation events
library(tidyverse)

#for working with CDS filtered
#variant.data<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/zm-sb_split_mafs/chr10.CDSonly.indelonly_reformated.tsv", col_names  = T)
Sb313.cds<-read_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb313.CDS.bed", col_names = c("CHROM","Start","End","ID","Strand"))

#for working with gene filtered
#variant.data<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/zm-sb_split_mafs/chr10.geneonly.indelonly_reformatted.tsv", col_names  = T)
#Sb313.gene<-read_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb313.gene.bed", col_names = c("CHROM","Start","End","ID","Strand"))

#bin1.variant.data<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/chr10.bin1.CDSonly.indelonly_reformatted.reducedheader.txt", col_names = T)
#bin2.variant.data<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/chr10.bin2.CDSonly.indelonly_reformatted.reducedheader.txt", col_names = T)

bin1.variant.data<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/test_mafs/chr10.bin1.CDSonly.indelonly_reformatted.reducedheader.txt", col_names = T)
bin2.variant.data<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/test_mafs/chr10.bin2.CDSonly.indelonly_reformatted.reducedheader.txt", col_names = T)


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

#split the alt allele to 1 and 2
#test to see it works
#head(variant.data) %>% 
#  mutate(ALT1 = str_split(string = ALT, pattern = ",", simplify = TRUE)[,1], 
#         ALT2 = str_split(string = ALT, pattern=",",simplify=TRUE)[,2]) %>% 
#  select(contains("ALT"))
#variant.data<-variant.data %>% mutate(ALT1 = str_split(string = ALT, pattern = ",", simplify = TRUE)[,1], 
#                        ALT2 = str_split(string = ALT, pattern=",",simplify=TRUE)[,2])
##where there isn't an ALT2, ALT2 == "" (not NA)

head(bin1.variant.data) %>% 
  mutate(ALT1 = str_split(string = ALT, pattern = ",", simplify = TRUE)[,1], 
                 ALT2 = str_split(string = ALT, pattern=",",simplify=TRUE)[,2]) %>% 
          select(contains("ALT"))
bin1.variant.data<-bin1.variant.data %>%  mutate(ALT1 = str_split(string = ALT, pattern = ",", simplify = TRUE)[,1], 
                                                 ALT2 = str_split(string = ALT, pattern=",",simplify=TRUE)[,2])
bin2.variant.data<-bin2.variant.data %>%  mutate(ALT1 = str_split(string = ALT, pattern = ",", simplify = TRUE)[,1], 
                                                 ALT2 = str_split(string = ALT, pattern=",",simplify=TRUE)[,2])
#IF THERE"S NOT AN ALT2 VARIANT, IT"S "" INSTEAD OF "NA"
filter(bin1.variant.data, ALT2 == "") %>% nrow() #(with the pre-CDS filtered vcf) 3161 out of 9611, just INDEL filter 76028/263859
filter(bin2.variant.data, ALT2 == "") %>% nrow() #(with the pre-CDS filtered vcf) 2396 out of 8221, just INDEL filter 39120/150811

#Filter out the insertions
#test to see it works
#head(variant.data, n=10) %>%
#  filter(str_length(REF) > str_length(ALT1) & str_length(REF) > str_length(ALT2))
#head(variant.data, n=10)

head(bin1.variant.data, n=10) %>% 
  filter(str_length(REF) > str_length(ALT1) | (str_length(REF) > str_length(ALT2) & ALT2 != ""))

#deletion.data<-filter(variant.data, str_length(REF) > str_length(ALT1) & str_length(REF) > str_length(ALT2))
#nrow(variant.data)-nrow(deletion.data)
#2870 insertions removed when using CDS filtered combined vcf
#15861 insertions removed when using gene filtered combined vcf

bin1.deletion.data<- filter(bin1.variant.data, str_length(REF) > str_length(ALT1) | (str_length(REF) > str_length(ALT2) & ALT2 != ""))
bin2.deletion.data<- filter(bin2.variant.data, str_length(REF) > str_length(ALT1) | (str_length(REF) > str_length(ALT2) & ALT2 != ""))
nrow(bin1.variant.data) - nrow(bin1.deletion.data) #1999 insertions removed when using CDS filtered bin1 vcf (50402 removed when using the INDEL only filtered vcf)
nrow(bin2.variant.data) - nrow(bin2.deletion.data) #1331 insertions removed when using CDS filtered bin2 vcf (25484 removed when using the INDEL only filtered vcf)

#Indicate which alleles are deletions and which are insertions
bin1.deletion.data<-mutate(bin1.deletion.data, ALT1_type = case_when(str_length(REF) > str_length(ALT1) & ALT1 != "*" ~ "deletion",
                                                 str_length(REF) < str_length(ALT1) ~ "insertion",
                                                 str_length(REF) == str_length(ALT1) ~ "not_indel",
                                                 ALT1 == "*" ~ "asterisk"),
       ALT2_type = case_when(str_length(REF) > str_length(ALT2) & ALT2 != "*" & ALT2 != "" ~ "deletion",
                             str_length(REF) < str_length(ALT2) ~ "insertion",
                             str_length(REF) == str_length(ALT2) ~ "not_indel",
                             ALT2 == "*" ~ "asterisk",
                             ALT2 == "" ~ "empty"))
bin2.deletion.data<-mutate(bin2.deletion.data, ALT1_type = case_when(str_length(REF) > str_length(ALT1) & ALT1 != "*" ~ "deletion",
                                                                     str_length(REF) < str_length(ALT1) ~ "insertion",
                                                                     str_length(REF) == str_length(ALT1) ~ "not_indel",
                                                                     ALT1 == "*" ~ "asterisk"),
                           ALT2_type = case_when(str_length(REF) > str_length(ALT2) & ALT2 != "*" & ALT2 != "" ~ "deletion",
                                                 str_length(REF) < str_length(ALT2) ~ "insertion",
                                                 str_length(REF) == str_length(ALT2) ~ "not_indel",
                                                 ALT2 == "*" ~ "asterisk",
                                                 ALT2 == "" ~ "empty"))


#Now I need to link deletions with CDS/gene IDs
#colnames(deletion.data)[1]<-"CHROM"
#colnames(Sb313.cds)
#colnames(Sb313.gene)
#deletion.data<-mutate(deletion.data, POSEnd = POS+1)
colnames(bin1.deletion.data)[1] <-"CHROM"
colnames(bin2.deletion.data)[1] <-"CHROM"

###I NEED TO MAKE SURE THAT THIS IS USING THE REF_SB313.CDS DF INSTEAD OF GENE###

#Overlap
#nrow(deletion.data) > nrow(Sb313.gene)
nrow(bin1.deletion.data) > nrow(ref_Sb313.cds) #want to get a false
nrow(bin2.deletion.data) > nrow(ref_Sb313.cds) #want to iterate through the shorter data frame

#deletion.data$SB_ID<-NA
bin1.deletion.data$CDS_ID <-NA
bin2.deletion.data$CDS_ID <-NA

#test
#deletion.data$SB_ID[1]<-filter(Sb313.cds, Start <= deletion.data$POS[1] & End >= deletion.data$POS[1]) %>% select(ID) %>% pull() %>% paste(collapse=":")
#deletion.data$SB_ID[1]<-filter(Sb313.gene, Start <= deletion.data$POS[1] & End >= deletion.data$POS[1]) %>% select(ID) %>% pull() %>% paste(collapse=":")

filter(ref_Sb313.cds, Start <= bin1.deletion.data$POS[11] & End >= bin1.deletion.data$POS[11] & CHROM == bin1.deletion.data$CHROM[11]) %>% select(CDS_ID) %>% pull()

#make sure to change to cds or gene depending on which filter you are using
#for(i in 1:nrow(deletion.data)){
#  deletion.data$SB_ID[i]<-filter(Sb313.gene, Start <= deletion.data$POS[i] & End >= deletion.data$POS[i]) %>% select(ID) %>% pull() %>% paste(collapse=":")
#}

for(i in 1:nrow(bin1.deletion.data)){
  t<-filter(ref_Sb313.cds, Start <= bin1.deletion.data$POS[i] & End >= bin1.deletion.data$POS[i] & CHROM == bin1.deletion.data$CHROM[i]) 
if(nrow(t) > 0){
  bin1.deletion.data$CDS_ID[i]<-select(t, CDS_ID) %>% pull()
}
}
bin1.deletion.data$CDS_ID %>% is.na() %>% summary() #see how many bin1 deletions were/weren't associated with a ref CDS
#1819 had a CDS ID, 5793 didn't (When using the CDS filtered vcf) 
#1819 had a CDS ID, 211638 didn't (when using the indel only filtered vcf) 
#[increased to 1897 vs. 218135 when not both alleles have to be deletion]
for(i in 1:nrow(bin2.deletion.data)){
  t<-filter(ref_Sb313.cds, Start <= bin2.deletion.data$POS[i] & End >= bin2.deletion.data$POS[i] & CHROM == bin2.deletion.data$CHROM[i]) 
  if(nrow(t) > 0){
    bin2.deletion.data$CDS_ID[i]<-select(t, CDS_ID) %>% pull()
  }
}
bin2.deletion.data$CDS_ID %>% is.na() %>% summary() #see how many bin2 deletions were/weren't associated with a ref CDS
#2037 had a CDS ID, 4853 didn't (when using the CDS filtered vcf)
#2037 had a CDS ID, 123290 didn't (when using the indel only filtered vcf)
#[increased to 2102 vs. 126485 when not both alleles have to be deletion]

#What are the lengths of these alt alleles (deletions)
bin1.deletion.data<-bin1.deletion.data %>% mutate(ALT1_Length = str_length(REF) - str_length(ALT1)) 
bin1.deletion.data<-bin1.deletion.data %>% mutate(ALT2_Length = str_length(REF) - str_length(ALT2)) 
bin2.deletion.data<-bin2.deletion.data %>% mutate(ALT1_Length = str_length(REF) - str_length(ALT1)) 
bin2.deletion.data<-bin2.deletion.data %>% mutate(ALT2_Length = str_length(REF) - str_length(ALT2)) 

for(i in 1:nrow(bin1.deletion.data)){
  if(bin1.deletion.data$ALT1[i] == "*"){
    bin1.deletion.data$ALT1_Length[i]<-"*"
  }
  if(bin1.deletion.data$ALT2[i] == "*"){
    bin1.deletion.data$ALT2_Length[i]<-"*"
  }
  if(bin1.deletion.data$ALT2[i] == ""){
    bin1.deletion.data$ALT2_Length[i] <- 0
  }
}

for(i in 1:nrow(bin2.deletion.data)){
  if(bin2.deletion.data$ALT1[i] == "*"){
    bin2.deletion.data$ALT1_Length[i]<-"*"
  }
  if(bin2.deletion.data$ALT2[i] == "*"){
    bin2.deletion.data$ALT2_Length[i]<-"*"
  }
  if(bin2.deletion.data$ALT2[i] == ""){
    bin2.deletion.data$ALT2_Length[i] <- 0
  }
}

bin1.ref.dels<-bin1.deletion.data %>% na.omit() #removes deletions without a CDS ID
bin2.ref.dels<-bin2.deletion.data %>% na.omit()

filter(bin1.ref.dels, ALT2 == "") %>% nrow() #619 out of 1820 don't have an ALT2 
filter(bin2.ref.dels, ALT2 == "") %>% nrow() #558 out of 2037 don't have an ALT2

ggplot(bin1.ref.dels, aes(x=ALT1_Length, y=ALT2_Length))+
  #geom_point(alpha = 0.5)+
  geom_jitter(alpha = 0.5)+
  scale_x_discrete(limits=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19",
                          "20","21","22","24","25","26","27","28","29","31","32","33","35","36","39","40","43","45",
                          "51","57","58","72","78","81","114","115","117","126","128","129","131","176","274","283","331",
                          "367","386","498","664","733","735","1176","1309","1419","1492","2368","3488","3491","3953","5631",
                          "5690","6700","6871","7581","12037","12476","19192","25593","26354","37999","42949","44854","*"))+
  scale_y_discrete(limits=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","20","21","22",
                            "23","24","26","27","29","31","33","36","40","43","126","127","206","456","732",'875',"909","1385",
                            "1415","1927","1943",'2432',"2382","2708","3145","3150","3319","4474","4808","4814","5792","7909",
                            "10199","12948","12977","18168","19168","38000","51419","*"))+
  theme(axis.text.x = element_text(angle = 90))

ggplot(bin2.ref.dels, aes(x=ALT1_Length, y=ALT2_Length))+
  #geom_point(alpha = 0.5)+
  geom_jitter(alpha = 0.5)+
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19",
                            "20","21","22","23","24","26","27","28","29","30","32","33","34","35","36","37","44","46","49",
                            "52","54","57","78","81","87","90", "98","129","156","193","221","228","278","282","290","358",
                            "423","476","498","513","633","640","679","686","732","802","832","935","1024","1061","1205","1210",
                            "1337","1410","1560","1894","2258","3114","4230","4251","4789","5683",
                            "6182","6601","7177","7885","8247","9012","11468","13381","13563","14178","19923","23167","25394",
                            "33204","38000","38228","44080","49902","*"))+
  scale_y_discrete(limits=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","21","22",
                            "23","24","27","30","31","39","42","45","51","59","78","81","96","173","179","190","202","234","363",
                            "407","498","508","560","620","820","1005","1846","1847","1886",
                            "1944",'2323',"2521","3481","4230","4242","4280","4647","4683","6136","8387","8421",
                            "11813","16402","23188","23505","23511","26834","29549","38000","41379","49798","51159","211881","*"))+
  theme(axis.text.x = element_text(angle = 90))

as.data.frame(table(bin1.ref.dels$CDS_ID)) %>%
  ggplot(aes(x=Freq))+
  geom_histogram(bins = 50)+
  ggtitle("Counts of deletion calls per CDS, Bin1")
as.data.frame(table(bin2.ref.dels$CDS_ID)) %>%
  ggplot(aes(x=Freq))+
  geom_histogram(bins = 50)+
  ggtitle("Counts of deletion calls per CDS, Bin2")

bin1.ref.dels<-inner_join(x=bin1.ref.dels, y=select(ref_Sb313.cds, c(CDS_ID,CDS_Length)), by="CDS_ID")
bin2.ref.dels<-inner_join(x=bin2.ref.dels, y=select(ref_Sb313.cds, c(CDS_ID,CDS_Length)), by="CDS_ID")

ggplot(bin1.ref.dels, aes(x=CDS_Length))+
  geom_point(aes(y=as.numeric(ALT2_Length)),color = "goldenrod",alpha=0.5, size=1)+
  geom_point(aes(y=as.numeric(ALT1_Length)),color = "darkblue",alpha=0.5, size=1)+
  geom_hline(yintercept=6511)+
  ylab("ALT Allele Length (deletion bp)")+ggtitle("Bin1 deletions")

ggplot(bin2.ref.dels, aes(x=CDS_Length))+
  geom_point(aes(y=as.numeric(ALT2_Length)),color = "goldenrod",alpha=0.5, size=1)+
  geom_point(aes(y=as.numeric(ALT1_Length)),color = "darkblue",alpha=0.5, size=1)+
  geom_hline(yintercept=6511)+
  ylab("ALT Allele Length (deletion bp)")+ggtitle("Bin2 deletions")

tibble(CDS_Length = c(bin1.ref.dels$CDS_Length,bin1.ref.dels$CDS_Length),
      ALT_Length = c(as.numeric(bin1.ref.dels$ALT1_Length), as.numeric(bin1.ref.dels$ALT2_Length))) %>%
  filter(ALT_Length <= 6511 & ALT_Length > 0) %>%
ggplot(aes(x=CDS_Length, y=ALT_Length))+
  geom_point(alpha=0.5, size=1)+
  ggtitle("Bin1 deletions")
tibble(CDS_Length = c(bin2.ref.dels$CDS_Length,bin2.ref.dels$CDS_Length),
       ALT_Length = c(as.numeric(bin2.ref.dels$ALT1_Length), as.numeric(bin2.ref.dels$ALT2_Length))) %>%
  filter(ALT_Length <= 6511 & ALT_Length > 0) %>%
  ggplot(aes(x=CDS_Length, y=ALT_Length))+
  geom_point(alpha=0.5, size=1)+
  ggtitle("Bin2 deletions")

#Colors to match Hufford et al 2012 Trends in Genetics Phylogeny colors
genome_colors<-c("Zm" = "goldenrod", 
                 "Zv" = "forestgreen",
                 "Zx" = "darkred",
                 "Zh" = "darkviolet",
                 "Zd" = "darkgrey",
                 "TdFL" = "darkorange")

###Assign fractionation status###
#if a CDS has a deletion of any kind, it should be considered fractionated
# 0 = ref (no deletion), 1 | 2 = alt (deletion), . = no call (NA)

#this makes a data frame where the genotypes are reduced to 0, 1, NA
#make insertion alleles and asterisk alleles NA
bin1.ref.dels.tmp<-mutate(bin1.ref.dels, TdFL.bin1 = case_when(Sb313_TdFL_Chr10.bin1 == "0" ~ 0,
                                            Sb313_TdFL_Chr10.bin1 == "1" & ALT1_type == "deletion" ~ 1,
                                            Sb313_TdFL_Chr10.bin1 == "2" & ALT2_type == "deletion" ~ 1,
                                            Sb313_TdFL_Chr10.bin1 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                                            Sb313_TdFL_Chr10.bin1 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                                            Sb313_TdFL_Chr10.bin1 == "." ~ NA),
       Zd.bin1 = case_when(Sb313_ZdGigi_Chr10.bin1 == "0" ~ 0,
                           Sb313_ZdGigi_Chr10.bin1 == "1" & ALT1_type == "deletion" ~ 1,
                           Sb313_ZdGigi_Chr10.bin1 == "2" & ALT2_type == "deletion" ~ 1,
                           Sb313_ZdGigi_Chr10.bin1 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                           Sb313_ZdGigi_Chr10.bin1 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                           Sb313_ZdGigi_Chr10.bin1 == "." ~ NA),
       Zv.bin1 = case_when(Sb313_ZvTIL01_Chr10.bin1 == "0" ~ 0,
                           Sb313_ZvTIL01_Chr10.bin1 == "1" & ALT1_type == "deletion" ~ 1,
                           Sb313_ZvTIL01_Chr10.bin1 == "2" & ALT2_type == "deletion" ~ 1,
                           Sb313_ZvTIL01_Chr10.bin1 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                           Sb313_ZvTIL01_Chr10.bin1 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                           Sb313_ZvTIL01_Chr10.bin1 == "." ~ NA),
       Zm.bin1 = case_when(Sb313_ZmB73_Chr10.bin1 == "0" ~ 0,
                           Sb313_ZmB73_Chr10.bin1 == "1" & ALT1_type == "deletion" ~ 1,
                           Sb313_ZmB73_Chr10.bin1 == "2" & ALT2_type == "deletion" ~ 1,
                           Sb313_ZmB73_Chr10.bin1 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                           Sb313_ZmB73_Chr10.bin1 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                           Sb313_ZmB73_Chr10.bin1 == "." ~ NA)) %>%
  select(c("CHROM","POS","CDS_ID","CDS_Length","TdFL.bin1","Zd.bin1","Zv.bin1","Zm.bin1"))

bin1.Uniq.CDS_ID<-bin1.ref.dels.tmp$CDS_ID %>% unique() #unique CDS IDs in bin1
#creating a coarse fractionation table for bin1 to collapse all deletions within a unique CDS ID together
#Also adds the variable Del_Cnts which is the number of deletions found in that particular CDS ID
bin1.coarse_fractionation<-tibble(CDS_ID = NA,
                                  CDS_Length = NA,
                                  Var_Cnts = NA,
                                  TdFL.bin1 = NA,
                                  Zd.bin1 = NA,
                                  Zv.bin1 = NA,
                                  Zm.bin1 = NA)
for(i in 1:length(bin1.Uniq.CDS_ID)){
  t<-filter(bin1.ref.dels.tmp, CDS_ID == bin1.Uniq.CDS_ID[i])
  c<-nrow(t) #number of deletions found in that CDS_ID
  bin1.coarse_fractionation<-add_row(bin1.coarse_fractionation, 
                                     CDS_ID = bin1.Uniq.CDS_ID[i],
                                     CDS_Length = t$CDS_Length %>% unique(),
                                     Var_Cnts = c,
                                     TdFL.bin1 = if(sum(t$TdFL.bin1 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$TdFL.bin1))){NA}else{0}
                                       },
                                     Zd.bin1 = if(sum(t$Zd.bin1 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$Zd.bin1))){NA}else{0}
                                     },
                                     Zv.bin1 = if(sum(t$Zv.bin1 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$Zv.bin1))){NA}else{0}
                                     },
                                     Zm.bin1 = if(sum(t$Zm.bin1 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$Zm.bin1))){NA}else{0}
                                     }
                                     )
}
bin1.coarse_fractionation<- bin1.coarse_fractionation[-1,]

#Redo the same process for bin2
#this makes a data frame where the genotypes are reduced to 0, 1, NA
bin2.ref.dels.tmp<-mutate(bin2.ref.dels, TdFL.bin2 = case_when(Sb313_TdFL_Chr10.bin2 == "0" ~ 0,
                                                               Sb313_TdFL_Chr10.bin2 == "1" & ALT1_type == "deletion" ~ 1,
                                                               Sb313_TdFL_Chr10.bin2 == "2" & ALT2_type == "deletion" ~ 1,
                                                               Sb313_TdFL_Chr10.bin2 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                                                               Sb313_TdFL_Chr10.bin2 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                                                               Sb313_TdFL_Chr10.bin2 == "." ~ NA),
                          Zd.bin2 = case_when(Sb313_ZdGigi_Chr10.bin2 == "0" ~ 0,
                                              Sb313_ZdGigi_Chr10.bin2 == "1" & ALT1_type == "deletion" ~ 1,
                                              Sb313_ZdGigi_Chr10.bin2 == "2" & ALT2_type == "deletion" ~ 1,
                                              Sb313_ZdGigi_Chr10.bin2 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                                              Sb313_ZdGigi_Chr10.bin2 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                                              Sb313_ZdGigi_Chr10.bin2 == "." ~ NA),
                          Zv.bin2 = case_when(Sb313_ZvTIL01_Chr10.bin2 == "0" ~ 0,
                                              Sb313_ZvTIL01_Chr10.bin2 == "1" & ALT1_type == "deletion" ~ 1,
                                              Sb313_ZvTIL01_Chr10.bin2 == "2" & ALT2_type == "deletion" ~ 1,
                                              Sb313_ZvTIL01_Chr10.bin2 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                                              Sb313_ZvTIL01_Chr10.bin2 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                                              Sb313_ZvTIL01_Chr10.bin2 == "." ~ NA),
                          Zm.bin2 = case_when(Sb313_ZmB73_Chr10.bin2 == "0" ~ 0,
                                              Sb313_ZmB73_Chr10.bin2 == "1" & ALT1_type == "deletion" ~ 1,
                                              Sb313_ZmB73_Chr10.bin2 == "2" & ALT2_type == "deletion" ~ 1,
                                              Sb313_ZmB73_Chr10.bin2 == "1" & (ALT1_type == "insertion" | ALT1_type == "asterisk") ~ NA, 
                                              Sb313_ZmB73_Chr10.bin2 == "2" & (ALT2_type == "insertion" | ALT2_type == "asterisk") ~ NA,
                                              Sb313_ZmB73_Chr10.bin2 == "." ~ NA)) %>%
  select(c("CHROM","POS","CDS_ID","CDS_Length","TdFL.bin2","Zd.bin2","Zv.bin2","Zm.bin2"))

bin2.Uniq.CDS_ID<-bin2.ref.dels.tmp$CDS_ID %>% unique() #unique CDS IDs in bin1
#creating a coarse fractionation table for bin1 to collapse all deletions within a unique CDS ID together
#Also adds the variable Del_Cnts which is the number of deletions found in that particular CDS ID
bin2.coarse_fractionation<-tibble(CDS_ID = NA,
                                  CDS_Length = NA,
                                  Var_Cnts = NA,
                                  TdFL.bin2 = NA,
                                  Zd.bin2 = NA,
                                  Zv.bin2 = NA,
                                  Zm.bin2 = NA)
for(i in 1:length(bin2.Uniq.CDS_ID)){
  t<-filter(bin2.ref.dels.tmp, CDS_ID == bin2.Uniq.CDS_ID[i])
  c<-nrow(t) #number of deletions found in that CDS_ID
  bin2.coarse_fractionation<-add_row(bin2.coarse_fractionation, 
                                     CDS_ID = bin2.Uniq.CDS_ID[i],
                                     CDS_Length = t$CDS_Length %>% unique(),
                                     Var_Cnts = c,
                                     TdFL.bin2 = if(sum(t$TdFL.bin2 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$TdFL.bin2))){NA}else{0}
                                     },
                                     Zd.bin2 = if(sum(t$Zd.bin2 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$Zd.bin2))){NA}else{0}
                                     },
                                     Zv.bin2 = if(sum(t$Zv.bin2 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$Zv.bin2))){NA}else{0}
                                     },
                                     Zm.bin2 = if(sum(t$Zm.bin2 > 0, na.rm = T)){1}else{
                                       if(all(is.na(t$Zm.bin2))){NA}else{0}
                                     }
  )
}
bin2.coarse_fractionation<- bin2.coarse_fractionation[-1,]

tibble(binID = c(rep("Bin1", nrow(bin1.coarse_fractionation)),rep("Bin2",nrow(bin2.coarse_fractionation))),
       Del_Cnts = c(bin1.coarse_fractionation$Del_Cnts, bin2.coarse_fractionation$Del_Cnts)) %>%
ggplot(aes(x=Del_Cnts))+
  geom_histogram(aes(fill = binID), alpha = 0.66, binwidth = 5)+
  ggtitle("Number of deletions called per CDS ID by Bin")

tibble(binID = c(rep("Bin1", nrow(bin1.coarse_fractionation)),rep("Bin2",nrow(bin2.coarse_fractionation))),
       Del_Cnts = c(bin1.coarse_fractionation$Del_Cnts, bin2.coarse_fractionation$Del_Cnts),
       CDS_Length = c(bin1.coarse_fractionation$CDS_Length, bin2.coarse_fractionation$CDS_Length)) %>%
  mutate(Del_Cnts_per_CDS_Length = Del_Cnts/CDS_Length) %>%
  ggplot(aes(x=Del_Cnts_per_CDS_Length))+
  geom_histogram(aes(fill = binID), alpha = 0.66)+
  ggtitle("Number of deletions / CDS Length per CDS ID")
  
coarse.fractionation<-full_join(x=bin1.coarse_fractionation, y=bin2.coarse_fractionation, by="CDS_ID", suffix = c(".bin1",".bin2"))
#write_tsv(coarse.fractionation, "/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/coarse.fractionation.table.tsv")

#originally had 867 rows (before relaxed del filtering)
#after relaxed del filtering (where ALT alleles didn't both have to be deletions), row count is 893
#if the CDS does not have deletions called in any genome for one bin, but does for the other
#the bin without deletions for that CDS will be assumed to have the ref allele and be retained across all genomes
#after the relaxed deletion filter, the number of CDS with NAs due to joins actually went up, not down...
###To test this assumption
bin2.nodel.CDS_ID<-filter(coarse.fractionation, is.na(CDS_Length.bin2)) %>% select(CDS_ID)
bin1.nodel.CDS_ID<-filter(coarse.fractionation, is.na(CDS_Length.bin1)) %>% select(CDS_ID)
filter(ref_Sb313.cds, CDS_ID %in% pull(bin2.nodel.CDS_ID)) %>% write_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/bin2.nodels.CDS.bed")
filter(ref_Sb313.cds, CDS_ID %in% pull(bin1.nodel.CDS_ID)) %>% write_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/bin1.nodels.CDS.bed")

###
#Using the above files, I can find the ones that 
#(1) weren't in the original GVCFs (likely misalignment, should be fractionated)
#(2) were in the original GVCFs, but didn't have a deletion called within those exons (should be retained)
###
TdFL.bin1.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin1.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
TdFL.bin1.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin1.retainedCDS.txt", col_names = "CDS_ID")
TdFL.bin2.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin2.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
TdFL.bin2.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin2.retainedCDS.txt", col_names = "CDS_ID")

Zd.bin1.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin1.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
Zd.bin1.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin1.retainedCDS.txt", col_names = "CDS_ID")
Zd.bin2.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin2.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
Zd.bin2.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin2.retainedCDS.txt", col_names = "CDS_ID")

Zv.bin1.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin1.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
Zv.bin1.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin1.retainedCDS.txt", col_names = "CDS_ID")
Zv.bin2.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin2.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
Zv.bin2.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin2.retainedCDS.txt", col_names = "CDS_ID")

Zm.bin1.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin1.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
Zm.bin1.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin1.retainedCDS.txt", col_names = "CDS_ID")
Zm.bin2.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin2.exonsNotInGVCF.txt", col_names = c("CDS_ID"))
Zm.bin2.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin2.retainedCDS.txt", col_names = "CDS_ID")

#Do it for TdFL
for(i in 1:nrow(coarse.fractionation)){ #for each row in coarse.fractionation
  if(coarse.fractionation$CDS_ID[i] %in% pull(TdFL.bin1.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"TdFL.bin1"])){ #if the CDS_ID is in list of exons not in GVCF and is NA
    coarse.fractionation[i,"TdFL.bin1"]<-1 #assign it as fractionated for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(TdFL.bin1.retained[,1]) & is.na(coarse.fractionation[i,"TdFL.bin1"])){ #if the CDS_ID is in list of exons in GVCF but without del
    coarse.fractionation[i,"TdFL.bin1"] <- 0 #assign it as retained for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(TdFL.bin2.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"TdFL.bin2"])){ #repeat the process for bin2
    coarse.fractionation[i,"TdFL.bin2"]<-1
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(TdFL.bin2.retained[,1]) & is.na(coarse.fractionation[i,"TdFL.bin2"])){
    coarse.fractionation[i,"TdFL.bin2"] <- 0
  }
}
#Zd
for(i in 1:nrow(coarse.fractionation)){ #for each row in coarse.fractionation
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zd.bin1.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"Zd.bin1"])){ #if the CDS_ID is in list of exons not in GVCF
    coarse.fractionation[i,"Zd.bin1"]<-1 #assign it as fractionated for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zd.bin1.retained[,1]) & is.na(coarse.fractionation[i,"Zd.bin1"])){ #if the CDS_ID is in list of exons in GVCF but without del
    coarse.fractionation[i,"Zd.bin1"] <- 0 #assign it as retained for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zd.bin2.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"Zd.bin2"])){ #repeat the process for bin2
    coarse.fractionation[i,"Zd.bin2"]<-1
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zd.bin2.retained[,1]) & is.na(coarse.fractionation[i,"Zd.bin2"])){
    coarse.fractionation[i,"Zd.bin2"] <- 0
  }
}

#Zv
for(i in 1:nrow(coarse.fractionation)){ #for each row in coarse.fractionation
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zv.bin1.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"Zv.bin1"])){ #if the CDS_ID is in list of exons not in GVCF
    coarse.fractionation[i,"Zv.bin1"]<-1 #assign it as fractionated for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zv.bin1.retained[,1]) & is.na(coarse.fractionation[i,"Zv.bin1"])){ #if the CDS_ID is in list of exons in GVCF but without del
    coarse.fractionation[i,"Zv.bin1"] <- 0 #assign it as retained for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zv.bin2.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"Zv.bin2"])){ #repeat the process for bin2
    coarse.fractionation[i,"Zv.bin2"]<-1
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zv.bin2.retained[,1]) & is.na(coarse.fractionation[i,"Zv.bin2"])){
    coarse.fractionation[i,"Zv.bin2"] <- 0
  }
}
#Zm
for(i in 1:nrow(coarse.fractionation)){ #for each row in coarse.fractionation
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zm.bin1.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"Zm.bin1"])){ #if the CDS_ID is in list of exons not in GVCF
    coarse.fractionation[i,"Zm.bin1"]<-1 #assign it as fractionated for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zm.bin1.retained[,1]) & is.na(coarse.fractionation[i,"Zm.bin1"])){ #if the CDS_ID is in list of exons in GVCF but without del
    coarse.fractionation[i,"Zm.bin1"] <- 0 #assign it as retained for that particular genome/homoeolog
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zm.bin2.exonsNotInGVCF[,1]) & is.na(coarse.fractionation[i,"Zm.bin2"])){ #repeat the process for bin2
    coarse.fractionation[i,"Zm.bin2"]<-1
  }
  if(coarse.fractionation$CDS_ID[i] %in% pull(Zm.bin2.retained[,1]) & is.na(coarse.fractionation[i,"Zm.bin2"])){
    coarse.fractionation[i,"Zm.bin2"] <- 0
  }
}

###

#How many of the reference CDS don't have any fractionation?
#Since this is preliminary and only looking at Chrom 10
##filter to only cds on chr10 then filter for those ids that aren't in the coarse.fractionation$CDS_ID
filter(ref_Sb313.cds, CHROM == "10") %>% filter(!CDS_ID %in% coarse.fractionation$CDS_ID) %>% nrow() #4915
#out of
filter(ref_Sb313.cds, CHROM == "10") %>% nrow() #5808

#How to add in those that are fully retained across both homoeologs in all genomes?
filter(ref_Sb313.cds, CHROM == "10") %>% filter(!CDS_ID %in% coarse.fractionation$CDS_ID) %>% write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/chr10.nodel.forEitherBin.bed",col_names = F)

##Add in those that weren't identified in any bin from the VCF

TdFL.bin1.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin1.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")
TdFL.bin2.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin2.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")

Zd.bin1.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin1.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")
Zd.bin2.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin2.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")

Zv.bin1.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin1.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")
Zv.bin2.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin2.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")

Zm.bin1.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin1.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")
Zm.bin2.extra.retained<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin2.nodel.forEitherBin_retainedCDS.txt", col_names = "CDS_ID")

extra.retained.df<-tibble(CDS_ID=TdFL.bin1.extra.retained$CDS_ID, TdFL.bin1 = 0)
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=TdFL.bin2.extra.retained$CDS_ID, TdFL.bin2 = 0), by = "CDS_ID")
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=Zd.bin1.extra.retained$CDS_ID, Zd.bin1 = 0), by = "CDS_ID")
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=Zd.bin2.extra.retained$CDS_ID, Zd.bin2 = 0), by = "CDS_ID")
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=Zv.bin1.extra.retained$CDS_ID, Zv.bin1 = 0), by = "CDS_ID")
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=Zv.bin2.extra.retained$CDS_ID, Zv.bin2 = 0), by = "CDS_ID")
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=Zm.bin1.extra.retained$CDS_ID, Zm.bin1 = 0), by = "CDS_ID")
extra.retained.df<-full_join(x=extra.retained.df, y=tibble(CDS_ID=Zm.bin2.extra.retained$CDS_ID, Zm.bin2 = 0), by = "CDS_ID")

#add in the extra retained
for(i in 1:nrow(extra.retained.df)){ #for each row in coarse.fractionation
  if(!extra.retained.df[i,1] %in% coarse.fractionation[,"CDS_ID"]){ #if the CDS from the list is not in the coarse fractionation df yet
    coarse.fractionation<-add_row(coarse.fractionation, CDS_ID = pull(extra.retained.df[i,1]),
                                  TdFL.bin1 = pull(extra.retained.df[i,"TdFL.bin1"]), TdFL.bin2 = pull(extra.retained.df[i,"TdFL.bin2"]),
                                  Zd.bin1 = pull(extra.retained.df[i,"Zd.bin1"]), Zd.bin2 = pull(extra.retained.df[i,"Zd.bin2"]),
                                  Zv.bin1 = pull(extra.retained.df[i,"Zv.bin1"]), Zv.bin2 = pull(extra.retained.df[i,"Zv.bin2"]),
                                  Zm.bin1 = pull(extra.retained.df[i,"Zm.bin1"]), Zm.bin2 = pull(extra.retained.df[i,"Zm.bin2"]),
                                  CDS_Length.bin1 = NA, Var_Cnts.bin1 = NA, CDS_Length.bin2=NA, Var_Cnts.bin2=NA) #make the new row with the information for TdFL.bin1, but NAs for everything else
  }
}

#How many of the reference CDS aren't included now that we've added in those we feel confident are retained?
#Since this is preliminary and only looking at Chrom 10
##filter to only cds on chr10 then filter for those ids that aren't in the coarse.fractionation$CDS_ID
filter(ref_Sb313.cds, CHROM == "10") %>% filter(!CDS_ID %in% coarse.fractionation$CDS_ID) %>% nrow() #818 aren't included
filter(ref_Sb313.cds, CHROM == "10") %>% filter(CDS_ID %in% coarse.fractionation$CDS_ID) %>% nrow() #4990 included 

#Let's try to add in the fractionated status for those exons that aren't in the GVCFs for both bins 
TdFL.bin1.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin1.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))
TdFL.bin2.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_TdFL_Chr10.bin2.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))

Zd.bin1.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin1.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))
Zd.bin2.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZdGigi_Chr10.bin2.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))

Zv.bin1.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin1.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))
Zv.bin2.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZvTIL01_Chr10.bin2.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))

Zm.bin1.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin1.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))
Zm.bin2.extra.exonsNotInGVCF<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/split_mafs_4genomeTrial/unzipped_gvcf/Sb313_ZmB73_Chr10.bin2.nodel.forEitherBin.exonsNotInGVCF.IDs.txt", col_names = c("CDS_ID"))

extra.exonsNotInGVCF.df<-tibble(CDS_ID=TdFL.bin1.extra.exonsNotInGVCF$CDS_ID, TdFL.bin1 = 1)
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=TdFL.bin2.extra.exonsNotInGVCF$CDS_ID, TdFL.bin2 = 1), by = "CDS_ID")
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=Zd.bin1.extra.exonsNotInGVCF$CDS_ID, Zd.bin1 = 1), by = "CDS_ID")
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=Zd.bin2.extra.exonsNotInGVCF$CDS_ID, Zd.bin2 = 1), by = "CDS_ID")
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=Zv.bin1.extra.exonsNotInGVCF$CDS_ID, Zv.bin1 = 1), by = "CDS_ID")
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=Zv.bin2.extra.exonsNotInGVCF$CDS_ID, Zv.bin2 = 1), by = "CDS_ID")
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=Zm.bin1.extra.exonsNotInGVCF$CDS_ID, Zm.bin1 = 1), by = "CDS_ID")
extra.exonsNotInGVCF.df<-full_join(x=extra.exonsNotInGVCF.df, y=tibble(CDS_ID=Zm.bin2.extra.exonsNotInGVCF$CDS_ID, Zm.bin2 = 1), by = "CDS_ID")

extra.exonsNotInGVCF.df<-mutate(extra.exonsNotInGVCF.df, TdFL.both = case_when(TdFL.bin1 == 1 & TdFL.bin2 == 1 ~ TRUE, .default = FALSE),
       Zd.both = case_when(Zd.bin1 == 1 & Zd.bin2 == 1 ~ TRUE, .default = FALSE),
       Zv.both = case_when(Zv.bin1 == 1 & Zv.bin2 == 1 ~ TRUE, .default = FALSE),
       Zm.both = case_when(Zm.bin1 == 1 & Zm.bin2 == 1 ~ TRUE, .default = FALSE))

for(i in 1:nrow(coarse.fractionation)){ #go row by row
  if(coarse.fractionation$CDS_ID[i] %in% extra.exonsNotInGVCF.df$CDS_ID){ #if the CDS ID is in the extra exons not in GVCF list
    df<-filter(extra.exonsNotInGVCF.df, CDS_ID == coarse.fractionation$CDS_ID[i]) #filter it out to just that exon
    if(is.na(coarse.fractionation$TdFL.bin1[i]) & is.na(coarse.fractionation$TdFL.bin2[i]) & df$TdFL.both == TRUE){ #if both bins are missing
      coarse.fractionation$TdFL.bin1[i]<-1
      coarse.fractionation$TdFL.bin2[i]<-1 # assign it as fractionated
    }
    if(is.na(coarse.fractionation$Zd.bin1[i]) & is.na(coarse.fractionation$Zd.bin2[i]) & df$Zd.both == TRUE){
      coarse.fractionation$Zd.bin1[i]<-1
      coarse.fractionation$Zd.bin2[i]<-1
    }
    if(is.na(coarse.fractionation$Zv.bin1[i]) & is.na(coarse.fractionation$Zv.bin2[i]) & df$Zv.both == TRUE){
      coarse.fractionation$Zv.bin1[i]<-1
      coarse.fractionation$Zv.bin2[i]<-1
    }
    if(is.na(coarse.fractionation$Zm.bin1[i]) & is.na(coarse.fractionation$Zm.bin2[i]) & df$Zm.both == TRUE){
      coarse.fractionation$Zm.bin1[i]<-1
      coarse.fractionation$Zm.bin2[i]<-1
    }
  }
}

#if the CDS_ID isn't in the coarse.fractionation already, we're going to ignore it for now




#RetentionStatus: Both Retained (0/0) | 
#Bin1 Retained Only (0/1) | Bin2 Retained Only (1/0) | 
#Bin1 Del/Bin2 Uncalled (1/NA) | Bin1 Uncalled/Bin2 (NA/1) | 
#Both lost (1/1) | Both uncalled (NA/NA)

coarse.fractionation<-mutate(coarse.fractionation, TdFLStatus = case_when(TdFL.bin1 == 0 & TdFL.bin2 == 0 ~ "Both_Retained",
                                                    TdFL.bin1 == 1 & TdFL.bin2 == 0 ~ "Bin2_Retained",
                                                    TdFL.bin1 == 0 & TdFL.bin2 == 1 ~ "Bin1_Retained",
                                                    TdFL.bin1 == 1 & TdFL.bin2 == 1 ~ "Both_Lost",
                                                    TdFL.bin1 == 1 & is.na(TdFL.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                                    TdFL.bin1 == 0 & is.na(TdFL.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                                    is.na(TdFL.bin1) & TdFL.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                                    is.na(TdFL.bin1) & TdFL.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                                    is.na(TdFL.bin1) & is.na(TdFL.bin2) ~ "Both_Uncalled"),
       ZdStatus = case_when(Zd.bin1 == 0 & Zd.bin2 == 0 ~ "Both_Retained",
                                Zd.bin1 == 1 & Zd.bin2 == 0 ~ "Bin2_Retained",
                                Zd.bin1 == 0 & Zd.bin2 == 1 ~ "Bin1_Retained",
                                Zd.bin1 == 1 & Zd.bin2 == 1 ~ "Both_Lost",
                                Zd.bin1 == 1 & is.na(Zd.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                Zd.bin1 == 0 & is.na(Zd.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                is.na(Zd.bin1) & Zd.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                is.na(Zd.bin1) & Zd.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                is.na(Zd.bin1) & is.na(Zd.bin2) ~ "Both_Uncalled"),
       ZvStatus = case_when(Zv.bin1 == 0 & Zv.bin2 == 0 ~ "Both_Retained",
                                 Zv.bin1 == 1 & Zv.bin2 == 0 ~ "Bin2_Retained",
                            Zv.bin1 == 0 & Zv.bin2 == 1 ~ "Bin1_Retained",
                            Zv.bin1 == 1 & Zv.bin2 == 1 ~ "Both_Lost",
                            Zv.bin1 == 1 & is.na(Zv.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                            Zv.bin1 == 0 & is.na(Zv.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                 is.na(Zv.bin1) & Zv.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                 is.na(Zv.bin1) & Zv.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                 is.na(Zv.bin1) & is.na(Zv.bin2) ~ "Both_Uncalled"),
       ZmStatus = case_when(Zm.bin1 == 0 & Zm.bin2 == 0 ~ "Both_Retained",
                            Zm.bin1 == 1 & Zm.bin2 == 0 ~ "Bin2_Retained",
                            Zm.bin1 == 0 & Zm.bin2 == 1 ~ "Bin1_Retained",
                            Zm.bin1 == 1 & Zm.bin2 == 1 ~ "Both_Lost",
                            Zm.bin1 == 1 & is.na(Zm.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                            Zm.bin1 == 0 & is.na(Zm.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                               is.na(Zm.bin1) & Zm.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                               is.na(Zm.bin1) & Zm.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                               is.na(Zm.bin1) & is.na(Zm.bin2) ~ "Both_Uncalled")
)

long.coarse.fractionation<-tibble(GeneID=c(rep(coarse.fractionation$CDS_ID, 4)),
                                        Genome=c(rep("TdFL",length(coarse.fractionation$CDS_ID)),
                                                 rep("Zd",length(coarse.fractionation$CDS_ID)),
                                                 rep("Zv",length(coarse.fractionation$CDS_ID)),
                                                 rep("Zm",length(coarse.fractionation$CDS_ID))),
                                        Status=c(coarse.fractionation$TdFLStatus,
                                                 coarse.fractionation$ZdStatus,
                                                 coarse.fractionation$ZvStatus,
                                                 coarse.fractionation$ZmStatus))
long.coarse.fractionation$Genome<-factor(long.coarse.fractionation$Genome, levels = c("TdFL","Zd","Zv","Zm"))
long.coarse.fractionation$Status<-factor(long.coarse.fractionation$Status, levels = c("Both_Retained", "Bin1_Retained","Bin2_Retained","Both_Lost",
                                                                                                  "Bin1_Retained:Bin2_Uncalled","Bin1_Lost:Bin2_Uncalled",
                                                                                                  "Bin1_Uncalled:Bin2_Retained","Bin1_Uncalled:Bin2_Lost",
                                                                                                  "Both_Uncalled"))
ggplot(data = long.coarse.fractionation, aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values = genome_colors)+
  xlab("Count By CDS")+
  ggtitle("Fractionation Status per unique CDS")

############################
#Let's think about the timing of deletion events
#For a deletion
#if it is shared by all genomes of bin1 | bin2 --> basal deletion
#if it is private to one genome of bin1 | bin2 --> species deletion
#if it is only in Zea genomes of bin1 | bin2 --> genus deletion
#if it is only in Zv and Zm of bin1 | bin2 --> sister deletion
#if it is any other sharing of bin1 | bin2 --> paraphyly deletion

coarse.fractionation<-mutate(coarse.fractionation, Bin1.Timing = case_when(
  TdFL.bin1 == 1 & Zd.bin1 == 1 & Zv.bin1 == 1 & Zm.bin1 == 1 ~ "basal",
  TdFL.bin1 == 0 & Zd.bin1 == 1 & Zv.bin1 == 1 & Zm.bin1 == 1 ~ "genus",
  TdFL.bin1 == 0 & Zd.bin1 == 0 & Zv.bin1 == 1 & Zm.bin1 == 1 ~ "sister",
  TdFL.bin1 == 1 & Zd.bin1 == 0 & Zv.bin1 == 0 & Zm.bin1 == 0 ~ "private",
  TdFL.bin1 == 0 & Zd.bin1 == 1 & Zv.bin1 == 0 & Zm.bin1 == 0 ~ "private",
  TdFL.bin1 == 0 & Zd.bin1 == 0 & Zv.bin1 == 1 & Zm.bin1 == 0 ~ "private",
  TdFL.bin1 == 0 & Zd.bin1 == 0 & Zv.bin1 == 0 & Zm.bin1 == 1 ~ "private",
  TdFL.bin1 == 1 & Zd.bin1 == 1 & Zv.bin1 == 1 & Zm.bin1 == 0 ~ "paraphyly",
  TdFL.bin1 == 1 & Zd.bin1 == 1 & Zv.bin1 == 0 & Zm.bin1 == 1 ~ "paraphyly",
  TdFL.bin1 == 1 & Zd.bin1 == 0 & Zv.bin1 == 1 & Zm.bin1 == 1 ~ "paraphyly",
  TdFL.bin1 == 1 & Zd.bin1 == 1 & Zv.bin1 == 0 & Zm.bin1 == 0 ~ "paraphyly",
  TdFL.bin1 == 1 & Zd.bin1 == 0 & Zv.bin1 == 1 & Zm.bin1 == 0 ~ "paraphyly",
  TdFL.bin1 == 1 & Zd.bin1 == 0 & Zv.bin1 == 0 & Zm.bin1 == 1 ~ "paraphyly",
  TdFL.bin1 == 0 & Zd.bin1 == 1 & Zv.bin1 == 1 & Zm.bin1 == 0 ~ "paraphyly",
  TdFL.bin1 == 0 & Zd.bin1 == 1 & Zv.bin1 == 0 & Zm.bin1 == 1 ~ "paraphyly",
  any(c(is.na(TdFL.bin1),is.na(Zd.bin1),is.na(Zv.bin1),is.na(Zm.bin1))) ~ "Undetermined"
),
Bin2.Timing = case_when(
  TdFL.bin2 == 1 & Zd.bin2 == 1 & Zv.bin2 == 1 & Zm.bin2 == 1 ~ "basal",
  TdFL.bin2 == 0 & Zd.bin2 == 1 & Zv.bin2 == 1 & Zm.bin2 == 1 ~ "genus",
  TdFL.bin2 == 0 & Zd.bin2 == 0 & Zv.bin2 == 1 & Zm.bin2 == 1 ~ "sister",
  TdFL.bin2 == 1 & Zd.bin2 == 0 & Zv.bin2 == 0 & Zm.bin2 == 0 ~ "private",
  TdFL.bin2 == 0 & Zd.bin2 == 1 & Zv.bin2 == 0 & Zm.bin2 == 0 ~ "private",
  TdFL.bin2 == 0 & Zd.bin2 == 0 & Zv.bin2 == 1 & Zm.bin2 == 0 ~ "private",
  TdFL.bin2 == 0 & Zd.bin2 == 0 & Zv.bin2 == 0 & Zm.bin2 == 1 ~ "private",
  TdFL.bin2 == 1 & Zd.bin2 == 1 & Zv.bin2 == 1 & Zm.bin2 == 0 ~ "paraphyly",
  TdFL.bin2 == 1 & Zd.bin2 == 1 & Zv.bin2 == 0 & Zm.bin2 == 1 ~ "paraphyly",
  TdFL.bin2 == 1 & Zd.bin2 == 0 & Zv.bin2 == 1 & Zm.bin2 == 1 ~ "paraphyly",
  TdFL.bin2 == 1 & Zd.bin2 == 1 & Zv.bin2 == 0 & Zm.bin2 == 0 ~ "paraphyly",
  TdFL.bin2 == 1 & Zd.bin2 == 0 & Zv.bin2 == 1 & Zm.bin2 == 0 ~ "paraphyly",
  TdFL.bin2 == 1 & Zd.bin2 == 0 & Zv.bin2 == 0 & Zm.bin2 == 1 ~ "paraphyly",
  TdFL.bin2 == 0 & Zd.bin2 == 1 & Zv.bin2 == 1 & Zm.bin2 == 0 ~ "paraphyly",
  TdFL.bin2 == 0 & Zd.bin2 == 1 & Zv.bin2 == 0 & Zm.bin2 == 1 ~ "paraphyly",
  any(c(is.na(TdFL.bin2),is.na(Zd.bin2),is.na(Zv.bin2),is.na(Zm.bin2))) ~ "Undetermined" 
  )
)
coarse.fractionation$Bin1.Timing<-factor(coarse.fractionation$Bin1.Timing, levels = c("basal","genus","sister","private","paraphyly","Undetermined"))
coarse.fractionation$Bin2.Timing<-factor(coarse.fractionation$Bin2.Timing, levels = c("basal","genus","sister","private","paraphyly","Undetermined"))

tibble(Timing = c(coarse.fractionation$Bin1.Timing, coarse.fractionation$Bin2.Timing),
       Bin = c(rep("Bin1",nrow(coarse.fractionation)),rep("Bin2",nrow(coarse.fractionation)))) %>%
  ggplot(aes(y=Timing))+
    geom_bar(aes(fill = Bin),position=position_dodge())+
    theme_minimal()+xlab("Count by CDS")+
    ggtitle("Timing using CDS")

#How much missing data is there for the Undetermined bins?





#Need to group by gene? 
deletion.data$GeneID<-NA
for(i in 1:nrow(deletion.data)){
  deletion.data$GeneID[i]<-deletion.data$SB_ID[i] %>% str_extract_all("ID.+?v3.1", simplify = T) %>% str_remove_all("ID=") %>% paste(collapse = ":")
}

deletion.data$GeneID %>% length() #19951 when using cds; 82266 when using gene
unique(deletion.data$GeneID) %>% length() #5428 when using cds; 9307 when using gene

#Let's add a deletionID (arbitrary) and then split GeneIDs so there's not multiple ID's per deletion
#then each gene would have a row and the deletion would be repeated

#this will take a long time to run (~2.5 hours)
deletion.data$Del_ID<-NA
for(i in 1:nrow(deletion.data)){
  deletion.data$Del_ID[i]<-paste("Del",i,sep="_")
}
deletionByGene<-tibble(CHROM=NA, POS=NA,ID=NA,REF=NA,
                       ALT=NA,QUAL=NA,FILTER=NA,INFO=NA,
                       FORMAT=NA,Sb313_TdFL_Chr10.bin1=NA,Sb313_TdFL_Chr10.bin2=NA,Sb313_ZdGigi_Chr10.bin1=NA,
                       Sb313_ZdGigi_Chr10.bin2=NA,Sb313_ZmB73_Chr10.bin1=NA,Sb313_ZmB73_Chr10.bin2=NA,Sb313_ZvTIL01_Chr10.bin1=NA,
                       Sb313_ZvTIL01_Chr10.bin2=NA,ALT1=NA,ALT2=NA,Del_ID=NA,GeneID=NA)
for(j in 1:nrow(deletion.data)){
  temp<-deletion.data[j,]
  for(i in str_split(temp$GeneID,pattern=":",simplify=TRUE)[1,]){
    deletionByGene[nrow(deletionByGene)+1,1:19]<-temp[1,1:19]
    deletionByGene[nrow(deletionByGene),20]<-temp[1,"Del_ID"]
    deletionByGene[nrow(deletionByGene),21]<-i
  }
  
}
deletionByGene<-deletionByGene[-1,]
#Let's make a coarse deletion by unique set of GeneIDs matrix
#Where if there's any deletion, the gene is considered fractionated
#Test
#temp<-filter(deletionByGene, GeneID == unique.GeneID[1])
#temp<-replace(temp,temp==".",NA)
#if(all(is.na(temp$Sb313_TdFL_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin2))) > 0){1}else{0}}

unique.GeneID <- deletionByGene$GeneID %>% unique()
coarse.fractionationByGeneID<-tibble(GeneID=NA,DelCnt=NA,TdFL.bin1=NA, TdFL.bin2=NA,
                                     ZdGigi.bin1=NA,ZdGigi.bin2=NA,
                                     ZmB73.bin1=NA,ZmB73.bin2=NA,
                                     ZvTIL01.bin1=NA,ZvTIL01.bin2=NA)
for(i in 1:length(unique.GeneID)){
  temp<-filter(deletionByGene, GeneID == unique.GeneID[i])
  temp<-replace(temp, temp==".",NA)
  TdFL.1<-if(all(is.na(temp$Sb313_TdFL_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin1))) > 0){1}else{0}}
  TdFL.2<-if(all(is.na(temp$Sb313_TdFL_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin2))) > 0){1}else{0}}
  ZdGigi.1<-if(all(is.na(temp$Sb313_ZdGigi_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZdGigi_Chr10.bin1))) > 0){1}else{0}}
  ZdGigi.2<-if(all(is.na(temp$Sb313_ZdGigi_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZdGigi_Chr10.bin2))) > 0){1}else{0}}
  ZmB73.1<-if(all(is.na(temp$Sb313_ZmB73_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZmB73_Chr10.bin1))) > 0){1}else{0}}
  ZmB73.2<-if(all(is.na(temp$Sb313_ZmB73_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZmB73_Chr10.bin2))) > 0){1}else{0}}
  ZvTIL01.1<-if(all(is.na(temp$Sb313_ZvTIL01_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZvTIL01_Chr10.bin1))) > 0){1}else{0}}
  ZvTIL01.2<-if(all(is.na(temp$Sb313_ZvTIL01_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZvTIL01_Chr10.bin2))) > 0){1}else{0}}
  coarse.fractionationByGeneID<-add_row(coarse.fractionationByGeneID, GeneID=unique.GeneID[i],DelCnt=length(unique(temp$Del_ID)),
                                        TdFL.bin1=TdFL.1, TdFL.bin2=TdFL.2,
                                        ZdGigi.bin1=ZdGigi.1,ZdGigi.bin2=ZdGigi.2,
                                        ZmB73.bin1=ZmB73.1,ZmB73.bin2=ZmB73.2,
                                        ZvTIL01.bin1=ZvTIL01.1,ZvTIL01.bin2=ZvTIL01.2)
}
coarse.fractionationByGeneID<-coarse.fractionationByGeneID[-1,] #removes starting line of NAs

ggplot(data = coarse.fractionationByGeneID, aes(x=DelCnt))+
  geom_histogram(binwidth = 10)+
  ggtitle("Unique Deletions per Ref. Gene ID")+
  xlab("Count of unique deletions")+ylab("")+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.DelsPerGene.png",device="png",dpi=300, height = 6, width = 6, units = "in")

mean(coarse.fractionationByGeneID$DelCnt, na.rm = T)
median(coarse.fractionationByGeneID$DelCnt, na.rm = T)
#Categorize Homoeolog "Fractionation" Status (COARSE)
#RetentionStatus: Both Retained (0/0) | Bin1 Retained Only (0/1) | Bin2 Retained Only (1/0) | Bin1 Del/Bin2 Uncalled (1/NA) | Bin1 Uncalled/Bin2 (NA/1) | Both lost (1/1) | Both uncalled (NA/NA)
coarse.fractionationByGeneID<-
  mutate(coarse.fractionationByGeneID, 
       TdFLStatus = case_when(TdFL.bin1 == 0 & TdFL.bin2 == 0 ~ "Both_Retained",
                              TdFL.bin1 == 1 & TdFL.bin2 == 0 ~ "Bin2_Retained",
                              TdFL.bin1 == 0 & TdFL.bin2 == 1 ~ "Bin1_Retained",
                              TdFL.bin1 == 1 & TdFL.bin2 == 1 ~ "Both_Lost",
                              TdFL.bin1 == 1 & is.na(TdFL.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                              TdFL.bin1 == 0 & is.na(TdFL.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                              is.na(TdFL.bin1) & TdFL.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                              is.na(TdFL.bin1) & TdFL.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                              is.na(TdFL.bin1) & is.na(TdFL.bin2) ~ "Both_Uncalled"),
       ZdGigiStatus = case_when(ZdGigi.bin1 == 0 & ZdGigi.bin2 == 0 ~ "Both_Retained",
                                ZdGigi.bin1 == 1 & ZdGigi.bin2 == 0 ~ "Bin2_Retained",
                                ZdGigi.bin1 == 0 & ZdGigi.bin2 == 1 ~ "Bin1_Retained",
                                ZdGigi.bin1 == 1 & ZdGigi.bin2 == 1 ~ "Both_Lost",
                                ZdGigi.bin1 == 1 & is.na(ZdGigi.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                ZdGigi.bin1 == 0 & is.na(ZdGigi.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                is.na(ZdGigi.bin1) & ZdGigi.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                is.na(ZdGigi.bin1) & ZdGigi.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                is.na(ZdGigi.bin1) & is.na(ZdGigi.bin2) ~ "Both_Uncalled"),
       ZvTIL01Status = case_when(ZvTIL01.bin1 == 0 & ZvTIL01.bin2 == 0 ~ "Both_Retained",
                                 ZvTIL01.bin1 == 1 & ZvTIL01.bin2 == 0 ~ "Bin2_Retained",
                                 ZvTIL01.bin1 == 0 & ZvTIL01.bin2 == 1 ~ "Bin1_Retained",
                                 ZvTIL01.bin1 == 1 & ZvTIL01.bin2 == 1 ~ "Both_Lost",
                                 ZvTIL01.bin1 == 1 & is.na(ZvTIL01.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                 ZvTIL01.bin1 == 0 & is.na(ZvTIL01.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                is.na(ZvTIL01.bin1) & ZvTIL01.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                is.na(ZvTIL01.bin1) & ZvTIL01.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                is.na(ZvTIL01.bin1) & is.na(ZvTIL01.bin2) ~ "Both_Uncalled"),
       ZmB73Status = case_when(ZmB73.bin1 == 0 & ZmB73.bin2 == 0 ~ "Both_Retained",
                               ZmB73.bin1 == 1 & ZmB73.bin2 == 0 ~ "Bin2_Retained",
                               ZmB73.bin1 == 0 & ZmB73.bin2 == 1 ~ "Bin1_Retained",
                               ZmB73.bin1 == 1 & ZmB73.bin2 == 1 ~ "Both_Lost",
                               ZmB73.bin1 == 1 & is.na(ZmB73.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                               ZmB73.bin1 == 0 & is.na(ZmB73.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                               is.na(ZmB73.bin1) & ZmB73.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                               is.na(ZmB73.bin1) & ZmB73.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                               is.na(ZmB73.bin1) & is.na(ZmB73.bin2) ~ "Both_Uncalled"),
       )
long.coarse.FractionationStatus<-tibble(GeneID=c(rep(coarse.fractionationByGeneID$GeneID, 4)),
                                        Genome=c(rep("TdFL",length(coarse.fractionationByGeneID$GeneID)),
                                                 rep("ZdGigi",length(coarse.fractionationByGeneID$GeneID)),
                                                 rep("ZvTIL01",length(coarse.fractionationByGeneID$GeneID)),
                                                 rep("ZmB73",length(coarse.fractionationByGeneID$GeneID))),
                                        Status=c(coarse.fractionationByGeneID$TdFLStatus,
                                                 coarse.fractionationByGeneID$ZdGigiStatus,
                                                 coarse.fractionationByGeneID$ZvTIL01Status,
                                                 coarse.fractionationByGeneID$ZmB73Status))
long.coarse.FractionationStatus$Genome<-factor(long.coarse.FractionationStatus$Genome, levels = c("TdFL","ZdGigi","ZvTIL01","ZmB73"))
long.coarse.FractionationStatus$Status<-factor(long.coarse.FractionationStatus$Status, levels = c("Both_Retained", "Bin1_Retained","Bin2_Retained","Both_Lost",
                                                                                                  "Bin1_Retained:Bin2_Uncalled","Bin1_Lost:Bin2_Uncalled",
                                                                                                  "Bin1_Uncalled:Bin2_Retained","Bin1_Uncalled:Bin2_Lost",
                                                                                                  "Both_Uncalled"))
ggplot(data = long.coarse.FractionationStatus, aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.CoarseFractStatus.png",device="png",dpi=300, height = 6, width = 6, units = "in")

#Both Retained doesn't mean much because this dataset by it's very nature has at least 1 deletion
filter(long.coarse.FractionationStatus, Status %in% c("Bin1_Retained","Bin2_Retained","Both_Lost")) %>%
  ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.CoarseFractStatusSubset.png",device="png",dpi=300, height = 6, width = 6, units = "in")

##### What about fractionation as size of deletion relative to gene length
deletionByGene[1,"REF"] %>% str_length()
deletionByGene<-mutate(deletionByGene, del_length = str_length(REF))

Sb313.gene<-mutate(Sb313.gene, 
       gene_length = End - Start,
       Gene_ID = str_extract_all(ID, "ID.+?v3.1", simplify = T) %>% str_remove_all("ID="))
unique(Sb313.gene$Gene_ID) %>% length()
nrow(Sb313.gene) 
unique.Sb313.gene<-unique(Sb313.gene$Gene_ID) 
unique.Sb313.gene.length<-tibble(Gene_ID=NA, Gene_length=NA)
for(i in 1:length(unique.Sb313.gene)){
  temp<-filter(Sb313.gene, Gene_ID == unique.Sb313.gene[i])
  unique.Sb313.gene.length<-add_row(unique.Sb313.gene.length, Gene_ID= unique.Sb313.gene[i], Gene_length = max(temp$gene_length))
}
unique.Sb313.gene.length<-unique.Sb313.gene.length[-1,]
#unique.Sb313.gene.length<- unique(unique.Sb313.gene.length)

deletionByGene<-left_join(x=deletionByGene, y=unique.Sb313.gene.length, by=c("GeneID" = "Gene_ID"))
deletionByGene<-mutate(deletionByGene, DelPropOfGene = del_length/Gene_length) 

filter(deletionByGene, DelPropOfGene <= 1.0) %>%
ggplot(aes(x=DelPropOfGene))+
  geom_histogram(binwidth = 0.01)
#Most deletions are very small relative to the size of the gene


filter(deletionByGene, DelPropOfGene >= 0.25) %>% nrow() / nrow(deletionByGene)
#look just at the deletions that are at least 25% the size of the gene length, that's 9587 unique deletions (3.96% of all deletions)
filter(deletionByGene, DelPropOfGene >= 0.25) %>% select("GeneID") %>% unique() %>% nrow()
#corresponding to 5227 unique genes

#look at fractionation status of the proportionally larger deletions
quarter.deletionByGene<-filter(deletionByGene, DelPropOfGene >= 0.25)
quarter.unique.GeneID <- filter(deletionByGene, DelPropOfGene >= 0.25) %>% select("GeneID") %>% unique() %>% pull()
quarter.coarse.fractionationByGeneID<-tibble(GeneID=NA,DelCnt=NA,TdFL.bin1=NA, TdFL.bin2=NA,
                                     ZdGigi.bin1=NA,ZdGigi.bin2=NA,
                                     ZmB73.bin1=NA,ZmB73.bin2=NA,
                                     ZvTIL01.bin1=NA,ZvTIL01.bin2=NA)

for(i in 1:length(quarter.unique.GeneID)){
  temp<-filter(quarter.deletionByGene, GeneID == quarter.unique.GeneID[i])
  temp<-replace(temp, temp==".",NA)
  TdFL.1<-if(all(is.na(temp$Sb313_TdFL_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin1))) > 0){1}else{0}}
  TdFL.2<-if(all(is.na(temp$Sb313_TdFL_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin2))) > 0){1}else{0}}
  ZdGigi.1<-if(all(is.na(temp$Sb313_ZdGigi_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZdGigi_Chr10.bin1))) > 0){1}else{0}}
  ZdGigi.2<-if(all(is.na(temp$Sb313_ZdGigi_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZdGigi_Chr10.bin2))) > 0){1}else{0}}
  ZmB73.1<-if(all(is.na(temp$Sb313_ZmB73_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZmB73_Chr10.bin1))) > 0){1}else{0}}
  ZmB73.2<-if(all(is.na(temp$Sb313_ZmB73_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZmB73_Chr10.bin2))) > 0){1}else{0}}
  ZvTIL01.1<-if(all(is.na(temp$Sb313_ZvTIL01_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZvTIL01_Chr10.bin1))) > 0){1}else{0}}
  ZvTIL01.2<-if(all(is.na(temp$Sb313_ZvTIL01_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZvTIL01_Chr10.bin2))) > 0){1}else{0}}
  quarter.coarse.fractionationByGeneID<-add_row(quarter.coarse.fractionationByGeneID, 
                                                GeneID=quarter.unique.GeneID[i],
                                                DelCnt=length(unique(temp$Del_ID)),
                                        TdFL.bin1=TdFL.1, TdFL.bin2=TdFL.2,
                                        ZdGigi.bin1=ZdGigi.1,ZdGigi.bin2=ZdGigi.2,
                                        ZmB73.bin1=ZmB73.1,ZmB73.bin2=ZmB73.2,
                                        ZvTIL01.bin1=ZvTIL01.1,ZvTIL01.bin2=ZvTIL01.2)
}
quarter.coarse.fractionationByGeneID<-quarter.coarse.fractionationByGeneID[-1,] #removes starting line of NAs

ggplot(data = quarter.coarse.fractionationByGeneID, aes(x=DelCnt))+
  geom_histogram(binwidth = 1)+
  ggtitle("Unique Deletions per Ref. Gene ID")+
  xlab("Count of unique deletions 25% gene length in size")+ylab("")+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.QuarterDelsPerGene.png",device="png",dpi=300, width = 6, height = 6, units = "in")

mean(quarter.coarse.fractionationByGeneID$DelCnt, na.rm = T)
median(quarter.coarse.fractionationByGeneID$DelCnt, na.rm = T)

quarter.coarse.fractionationByGeneID<-
  mutate(quarter.coarse.fractionationByGeneID, 
         TdFLStatus = case_when(TdFL.bin1 == 0 & TdFL.bin2 == 0 ~ "Both_Retained",
                                TdFL.bin1 == 1 & TdFL.bin2 == 0 ~ "Bin2_Retained",
                                TdFL.bin1 == 0 & TdFL.bin2 == 1 ~ "Bin1_Retained",
                                TdFL.bin1 == 1 & TdFL.bin2 == 1 ~ "Both_Lost",
                                TdFL.bin1 == 1 & is.na(TdFL.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                TdFL.bin1 == 0 & is.na(TdFL.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                is.na(TdFL.bin1) & TdFL.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                is.na(TdFL.bin1) & TdFL.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                is.na(TdFL.bin1) & is.na(TdFL.bin2) ~ "Both_Uncalled"),
         ZdGigiStatus = case_when(ZdGigi.bin1 == 0 & ZdGigi.bin2 == 0 ~ "Both_Retained",
                                  ZdGigi.bin1 == 1 & ZdGigi.bin2 == 0 ~ "Bin2_Retained",
                                  ZdGigi.bin1 == 0 & ZdGigi.bin2 == 1 ~ "Bin1_Retained",
                                  ZdGigi.bin1 == 1 & ZdGigi.bin2 == 1 ~ "Both_Lost",
                                  ZdGigi.bin1 == 1 & is.na(ZdGigi.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                  ZdGigi.bin1 == 0 & is.na(ZdGigi.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                  is.na(ZdGigi.bin1) & ZdGigi.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                  is.na(ZdGigi.bin1) & ZdGigi.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                  is.na(ZdGigi.bin1) & is.na(ZdGigi.bin2) ~ "Both_Uncalled"),
         ZvTIL01Status = case_when(ZvTIL01.bin1 == 0 & ZvTIL01.bin2 == 0 ~ "Both_Retained",
                                   ZvTIL01.bin1 == 1 & ZvTIL01.bin2 == 0 ~ "Bin2_Retained",
                                   ZvTIL01.bin1 == 0 & ZvTIL01.bin2 == 1 ~ "Bin1_Retained",
                                   ZvTIL01.bin1 == 1 & ZvTIL01.bin2 == 1 ~ "Both_Lost",
                                   ZvTIL01.bin1 == 1 & is.na(ZvTIL01.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                   ZvTIL01.bin1 == 0 & is.na(ZvTIL01.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                   is.na(ZvTIL01.bin1) & ZvTIL01.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                   is.na(ZvTIL01.bin1) & ZvTIL01.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                   is.na(ZvTIL01.bin1) & is.na(ZvTIL01.bin2) ~ "Both_Uncalled"),
         ZmB73Status = case_when(ZmB73.bin1 == 0 & ZmB73.bin2 == 0 ~ "Both_Retained",
                                 ZmB73.bin1 == 1 & ZmB73.bin2 == 0 ~ "Bin2_Retained",
                                 ZmB73.bin1 == 0 & ZmB73.bin2 == 1 ~ "Bin1_Retained",
                                 ZmB73.bin1 == 1 & ZmB73.bin2 == 1 ~ "Both_Lost",
                                 ZmB73.bin1 == 1 & is.na(ZmB73.bin2) ~ "Bin1_Lost:Bin2_Uncalled",
                                 ZmB73.bin1 == 0 & is.na(ZmB73.bin2) ~ "Bin1_Retained:Bin2_Uncalled",
                                 is.na(ZmB73.bin1) & ZmB73.bin2 == 1 ~ "Bin1_Uncalled:Bin2_Lost",
                                 is.na(ZmB73.bin1) & ZmB73.bin2 == 0 ~ "Bin1_Uncalled:Bin2_Retained",
                                 is.na(ZmB73.bin1) & is.na(ZmB73.bin2) ~ "Both_Uncalled"),
  )
long.quarter.coarse.FractionationStatus<-tibble(GeneID=c(rep(quarter.coarse.fractionationByGeneID$GeneID, 4)),
                                        Genome=c(rep("TdFL",length(quarter.coarse.fractionationByGeneID$GeneID)),
                                                 rep("ZdGigi",length(quarter.coarse.fractionationByGeneID$GeneID)),
                                                 rep("ZvTIL01",length(quarter.coarse.fractionationByGeneID$GeneID)),
                                                 rep("ZmB73",length(quarter.coarse.fractionationByGeneID$GeneID))),
                                        Status=c(quarter.coarse.fractionationByGeneID$TdFLStatus,
                                                 quarter.coarse.fractionationByGeneID$ZdGigiStatus,
                                                 quarter.coarse.fractionationByGeneID$ZvTIL01Status,
                                                 quarter.coarse.fractionationByGeneID$ZmB73Status))
long.quarter.coarse.FractionationStatus$Genome<-factor(long.quarter.coarse.FractionationStatus$Genome, levels = c("TdFL","ZdGigi","ZvTIL01","ZmB73"))
long.quarter.coarse.FractionationStatus$Status<-factor(long.quarter.coarse.FractionationStatus$Status, levels = c("Both_Retained", "Bin1_Retained","Bin2_Retained","Both_Lost",
                                                                                                  "Bin1_Retained:Bin2_Uncalled","Bin1_Lost:Bin2_Uncalled",
                                                                                                  "Bin1_Uncalled:Bin2_Retained","Bin1_Uncalled:Bin2_Lost",
                                                                                                  "Both_Uncalled"))

ggplot(data = long.quarter.coarse.FractionationStatus, aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.QuarterCoarseFractStatus.png",device="png",dpi=300, width = 6, height=6, units = "in")

#Both Retained doesn't mean much because this dataset by it's very nature has at least 1 deletion
filter(long.quarter.coarse.FractionationStatus, Status %in% c("Both_Retained","Bin1_Retained","Bin2_Retained","Both_Lost")) %>%
  ggplot(aes(y=Status))+
  geom_bar(aes(fill=Genome), position = position_dodge())+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.Quarter.CoarseFractStatusSubset.png",device="png",dpi=300, width = 6, height = 6, units = "in")


############################
#Let's think about the timing of deletion events
#For a deletion
#if it is shared by all genomes of bin1 | bin2 --> basal deletion
#if it is private to one genome of bin1 | bin2 --> species deletion
#if it is only in Zea genomes of bin1 | bin2 --> genus deletion
#if it is only in Zv and Zm of bin1 | bin2 --> sister deletion
#if it is any other sharing of bin1 | bin2 --> paraphyly deletion

deletion.data<-mutate(deletion.data, Del_length = str_length(REF))

del.timing<-tibble(Del_ID = NA, Del_length = NA, GeneID = NA, 
                   TdFL.bin1=NA,TdFL.bin2=NA,ZdGigi.bin1=NA,ZdGigi.bin2=NA,
                   ZvTIL01.bin1=NA,ZvTIL01.bin2=NA,ZmB73.bin1=NA,ZmB73.bin2=NA)
for(i in 1:nrow(deletion.data)){
  temp<-deletion.data[i,]
  temp<-replace(temp, temp==".",NA)
  TdFL.1<-if(all(is.na(temp$Sb313_TdFL_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin1))) > 0){1}else{0}}
  TdFL.2<-if(all(is.na(temp$Sb313_TdFL_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_TdFL_Chr10.bin2))) > 0){1}else{0}}
  ZdGigi.1<-if(all(is.na(temp$Sb313_ZdGigi_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZdGigi_Chr10.bin1))) > 0){1}else{0}}
  ZdGigi.2<-if(all(is.na(temp$Sb313_ZdGigi_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZdGigi_Chr10.bin2))) > 0){1}else{0}}
  ZmB73.1<-if(all(is.na(temp$Sb313_ZmB73_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZmB73_Chr10.bin1))) > 0){1}else{0}}
  ZmB73.2<-if(all(is.na(temp$Sb313_ZmB73_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZmB73_Chr10.bin2))) > 0){1}else{0}}
  ZvTIL01.1<-if(all(is.na(temp$Sb313_ZvTIL01_Chr10.bin1))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZvTIL01_Chr10.bin1))) > 0){1}else{0}}
  ZvTIL01.2<-if(all(is.na(temp$Sb313_ZvTIL01_Chr10.bin2))){NA}else{if(sum(na.omit(as.numeric(temp$Sb313_ZvTIL01_Chr10.bin2))) > 0){1}else{0}}
  del.timing<-add_row(del.timing,Del_ID=temp$Del_ID[1],Del_length=temp$Del_length[1],GeneID=temp$GeneID[1],
                      TdFL.bin1=TdFL.1,TdFL.bin2=TdFL.2,ZdGigi.bin1=ZdGigi.1,ZdGigi.bin2=ZdGigi.2,
                      ZvTIL01.bin1=ZvTIL01.1,ZvTIL01.bin2=ZvTIL01.2,ZmB73.bin1=ZmB73.1,ZmB73.bin2=ZmB73.2)
}
del.timing<-del.timing[-1,]

del.timing<-mutate(del.timing, 
       Del_time.bin1 = case_when(TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 1 ~ "basal",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 0 ~ "private",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 0 ~ "private",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "private",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "private",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 1 ~ "genus",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 1 ~ "sister",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "paraphyly",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "paraphyly",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 1 ~ "paraphyly",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 0 ~ "paraphyly",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "paraphyly",
                                 TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "paraphyly",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "paraphyly",
                                 TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "paraphyly"),
       Del_time.bin2 = case_when(TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 1 ~ "basal",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 0 ~ "private",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 0 ~ "private",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "private",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "private",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 1 ~ "genus",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 1 ~ "sister",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "paraphyly",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "paraphyly",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 1 ~ "paraphyly",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 0 ~ "paraphyly",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "paraphyly",
                                 TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "paraphyly",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "paraphyly",
                                 TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "paraphyly"
       )) 

filter(del.timing, !is.na(Del_time.bin1) & !is.na(Del_time.bin2)) #so there are some deletions that are across both bins? Yes, that would be the multiple alt alleles at play

long.del.timing<-tibble(Del_ID = rep(del.timing$Del_ID, 2),
                        Del_length = rep(del.timing$Del_length, 2),
                        GeneID = rep(del.timing$GeneID, 2),
                        Bin = c(rep("bin1", nrow(del.timing)),rep("bin2",nrow(del.timing))),
                        Timing = c(del.timing$Del_time.bin1, del.timing$Del_time.bin2)
                        )
long.del.timing$Timing <- factor(long.del.timing$Timing, levels = c("basal","genus","sister","private","paraphyly"))

ggplot(long.del.timing, aes(y=Timing))+
  geom_bar(aes(fill=Bin),position=position_dodge())+
  theme_minimal()+
  ggtitle("Sharing of deletions across genomes")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.CoarseTiming.png",device="png",dpi=300, width = 6, height = 6, units="in")

filter(long.del.timing, Del_length > 100) %>%
  ggplot(aes(y=Timing))+
  geom_bar(aes(fill=Bin),position=position_dodge())+
  theme_minimal()+
  ggtitle("Sharing of deletions > 100bp across genomes")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.100bpDels.CoarseTiming.png",device="png",dpi=300, width = 6, height = 6, units="in")

long.del.timing %>% count(Bin, Timing)
filter(long.del.timing, Del_length > 100) %>% count(Bin, Timing)

paraphyly.dels<-filter(del.timing, Del_time.bin1 == "paraphyly" | Del_time.bin2 == "paraphyly")
paraphyly.dels<- mutate(paraphyly.dels, 
       bin1.paraphyly = case_when(TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "Td-Zd-Zv",
                                  TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "Td-Zd-Zm",
                                  TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 1 ~ "Td-Zv-Zm",
                                  TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 0 ~ "Td-Zd",
                                  TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "Td-Zv",
                                  TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "Td-Zm",
                                  TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "Zd-Zv",
                                  TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "Zd-Zm"),
       bin2.paraphyly = case_when(TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "Td-Zd-Zv",
                                  TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "Td-Zd-Zm",
                                  TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 1 ~ "Td-Zv-Zm",
                                  TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 0 ~ "Td-Zd",
                                  TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "Td-Zv",
                                  TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "Td-Zm",
                                  TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "Zd-Zv",
                                  TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "Zd-Zm")
       ) 
long.paraphyly.dels<-tibble(Del_ID = rep(paraphyly.dels$Del_ID, 2),
                            Del_length = rep(paraphyly.dels$Del_length, 2),
                            GeneID = rep(paraphyly.dels$GeneID, 2),
                            Bin = c(rep("bin1", nrow(paraphyly.dels)),rep("bin2",nrow(paraphyly.dels))),
                            Paraphyly = c(paraphyly.dels$bin1.paraphyly, paraphyly.dels$bin2.paraphyly)
)
long.paraphyly.dels<- na.omit(long.paraphyly.dels)
long.paraphyly.dels$Paraphyly<-factor(long.paraphyly.dels$Paraphyly, levels = c("Td-Zd-Zv","Td-Zd-Zm","Td-Zv-Zm","Td-Zd","Td-Zv","Td-Zm","Zd-Zv","Zd-Zm"))
ggplot(long.paraphyly.dels, aes(y=Paraphyly))+
  geom_bar(aes(fill=Bin), position=position_dodge())+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.TypesOfParaphyly.png",device="png",dpi=300, width = 6, height = 6, units="in")

mean(deletion.data$Del_length, na.rm = T)
median(deletion.data$Del_length, na.rm = T)

long.paraphyly.dels %>% count(Bin, Paraphyly)

paraphyly.100bp.dels<-filter(del.timing, Del_length > 100) %>% filter(Del_time.bin1 == "paraphyly" | Del_time.bin2 == "paraphyly")
paraphyly.100bp.dels<- mutate(paraphyly.100bp.dels, 
                        bin1.paraphyly = case_when(TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "Td-Zd-Zv",
                                                   TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "Td-Zd-Zm",
                                                   TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 1 ~ "Td-Zv-Zm",
                                                   TdFL.bin1 == 1 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 0 ~ "Td-Zd",
                                                   TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "Td-Zv",
                                                   TdFL.bin1 == 1 & ZdGigi.bin1 == 0 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "Td-Zm",
                                                   TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 1 & ZmB73.bin1 == 0 ~ "Zd-Zv",
                                                   TdFL.bin1 == 0 & ZdGigi.bin1 == 1 & ZvTIL01.bin1 == 0 & ZmB73.bin1 == 1 ~ "Zd-Zm"),
                        bin2.paraphyly = case_when(TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "Td-Zd-Zv",
                                                   TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "Td-Zd-Zm",
                                                   TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 1 ~ "Td-Zv-Zm",
                                                   TdFL.bin2 == 1 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 0 ~ "Td-Zd",
                                                   TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "Td-Zv",
                                                   TdFL.bin2 == 1 & ZdGigi.bin2 == 0 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "Td-Zm",
                                                   TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 1 & ZmB73.bin2 == 0 ~ "Zd-Zv",
                                                   TdFL.bin2 == 0 & ZdGigi.bin2 == 1 & ZvTIL01.bin2 == 0 & ZmB73.bin2 == 1 ~ "Zd-Zm")
) 
long.paraphyly.100bp.dels<-tibble(Del_ID = rep(paraphyly.100bp.dels$Del_ID, 2),
                            Del_length = rep(paraphyly.100bp.dels$Del_length, 2),
                            GeneID = rep(paraphyly.100bp.dels$GeneID, 2),
                            Bin = c(rep("bin1", nrow(paraphyly.100bp.dels)),rep("bin2",nrow(paraphyly.100bp.dels))),
                            Paraphyly = c(paraphyly.100bp.dels$bin1.paraphyly, paraphyly.100bp.dels$bin2.paraphyly)
)
long.paraphyly.100bp.dels<- na.omit(long.paraphyly.100bp.dels)
long.paraphyly.100bp.dels$Paraphyly<-factor(long.paraphyly.100bp.dels$Paraphyly, levels = c("Td-Zd-Zv","Td-Zd-Zm","Td-Zv-Zm","Td-Zd","Td-Zv","Td-Zm","Zd-Zv","Zd-Zm"))

ggplot(long.paraphyly.100bp.dels, aes(y=Paraphyly))+
  geom_bar(aes(fill=Bin), position=position_dodge())+
  theme_minimal()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Zeavolution.TypesOfParaphyly.100bp.png",device="png",dpi=300, width = 6, height = 6, units="in")

long.paraphyly.100bp.dels %>% count(Bin, Paraphyly)
