####Libraries to load####
library(tidyverse)

####Handy vectors####

#tripsacinae_genomes_IDs<-c("TdKS","TdFL","ZnPI615697","ZdGigi","ZdMomo","ZnPI615697_4to1","ZdGigi_4to1","ZdMomo_4to1","ZhRIMHU001","ZxTIL18","ZxTIL25","ZvTIL01","ZvTIL11",...)
#genome_colors
#chr_colors

####inputs to load####

#Sorghum CDS references
#ref_Sb313.cds<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed", col_names = c("CHROM","Start","End","ID","Strand","CDS_ID","Gene_ID","CDS_Length"))

#Headerless VCFs per Sb ref chr per subgenome

for(i in 1:10){
  for(b in c("bin1", "bin2")){
    if(i != 10){
      assign(paste0("chr0",i,"_",b,"_collapseddels"), 
             read_tsv(file = paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/noZn_vcf/chr0",i,".",b,".del.exonic.collapsed.bed"), 
                      col_names = c("CHROM","Start","Stop","REF","ALT","QUAL","TdFL","TdKS","ZdGigi_4to1",
                                   "ZdMomo_4to1","ZhRIMHU001","ZmB73","ZmB97","ZmCML103","ZmCML228",
                                   "ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmCML52","ZmCML69",
                                   "ZmHP301", "ZmIL14H", "ZmKi11","ZmKi3","ZmKy21", "ZmM162W","ZmM37W",
                                   "ZmMS71","ZmMo18W","ZmNC350","ZmNC358","ZmOh43","ZmOh7b",
                                   "ZmP39","ZmTx303","ZmTzi8","ZvTIL01","ZvTIL11","ZxTIL18","ZxTIL25"),
                      na = c("",NA,"."))
      )
    }else{
      assign(paste0("chr",i,"_",b,"_collapseddels"), 
             read_tsv(file = paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/noZn_vcf/chr",i,".",b,".del.exonic.collapsed.bed"), 
                      col_names = c("CHROM","Start","Stop","REF","ALT","QUAL","TdFL","TdKS","ZdGigi_4to1",
                                   "ZdMomo_4to1","ZhRIMHU001","ZmB73","ZmB97","ZmCML103","ZmCML228",
                                   "ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmCML52","ZmCML69",
                                   "ZmHP301", "ZmIL14H", "ZmKi11","ZmKi3","ZmKy21", "ZmM162W","ZmM37W",
                                   "ZmMS71","ZmMo18W","ZmNC350","ZmNC358","ZmOh43","ZmOh7b",
                                   "ZmP39","ZmTx303","ZmTzi8","ZvTIL01","ZvTIL11","ZxTIL18","ZxTIL25"),
                      na = c("",NA,".")) 
      )
    }
  }
}
####Combine collapsed deletion objects into one object####

ALL_collapseddels<-tibble(CHROM=NA,Start=NA,Stop=NA,REF=NA,ALT=NA,QUAL=NA,M=NA,TdFL=NA,TdKS=NA,ZdGigi_4to1=NA,
                          ZdMomo_4to1=NA,ZhRIMHU001=NA,ZmB73=NA,ZmB97=NA,ZmCML103=NA,ZmCML228=NA,
                          ZmCML247=NA,ZmCML277=NA,ZmCML322=NA,ZmCML333=NA,ZmCML52=NA,ZmCML69=NA,
                          ZmHP301=NA,ZmIL14H=NA,ZmKi11=NA,ZmKi3=NA,ZmKy21=NA,ZmM162W=NA,ZmM37W=NA,
                          ZmMS71=NA,ZmMo18W=NA,ZmNC350=NA,ZmNC358=NA,ZmOh43=NA,ZmOh7b=NA,
                          ZmP39=NA,ZmTx303=NA,ZmTzi8=NA,ZvTIL01=NA,ZvTIL11=NA,ZxTIL18=NA,ZxTIL25=NA)
for(i in ls(pattern = "collapseddels")){
  if(str_detect(i,"bin1")){
    ALL_collapseddels<-add_row(ALL_collapseddels,
                               CHROM=as.numeric(pull(get(i)[,"CHROM"])),Start=pull(get(i)[,"Start"]),Stop=pull(get(i)[,"Stop"]),
                               REF=pull(get(i)[,"REF"]),ALT=pull(get(i)[,"ALT"]),QUAL=pull(get(i)[,"QUAL"]),M="M1",
                               TdFL=pull(get(i)[,"TdFL"]),TdKS=pull(get(i)[,"TdKS"]),ZdGigi_4to1=pull(get(i)[,"ZdGigi_4to1"]),
                               ZdMomo_4to1=pull(get(i)[,"ZdMomo_4to1"]),ZhRIMHU001=pull(get(i)[,"ZhRIMHU001"]),
                               ZmB73=pull(get(i)[,"ZmB73"]),
                               ZmB97=pull(get(i)[,"ZmB97"]),ZmCML103=pull(get(i)[,"ZmCML103"]),ZmCML228=pull(get(i)[,"ZmCML228"]),
                               ZmCML247=pull(get(i)[,"ZmCML247"]),ZmCML277=pull(get(i)[,"ZmCML277"]),ZmCML322=pull(get(i)[,"ZmCML322"]),
                               ZmCML333=pull(get(i)[,"ZmCML333"]),ZmCML52=pull(get(i)[,"ZmCML52"]),ZmCML69=pull(get(i)[,"ZmCML69"]),
                               ZmHP301=pull(get(i)[,"ZmHP301"]),ZmIL14H=pull(get(i)[,"ZmIL14H"]),ZmKi11=pull(get(i)[,"ZmKi11"]),
                               ZmKi3=pull(get(i)[,"ZmKi3"]),ZmKy21=pull(get(i)[,"ZmKy21"]),ZmM162W=pull(get(i)[,"ZmM162W"]),
                               ZmM37W=pull(get(i)[,"ZmM37W"]), ZmMS71=pull(get(i)[,"ZmMS71"]),ZmMo18W=pull(get(i)[,"ZmMo18W"]),
                               ZmNC350=pull(get(i)[,"ZmNC350"]),ZmNC358=pull(get(i)[,"ZmNC358"]),ZmOh43=pull(get(i)[,"ZmOh43"]),
                               ZmOh7b=pull(get(i)[,"ZmOh7b"]), ZmP39=pull(get(i)[,"ZmP39"]),ZmTx303=pull(get(i)[,"ZmTx303"]),
                               ZmTzi8=pull(get(i)[,"ZmTzi8"]),ZvTIL01=pull(get(i)[,"ZvTIL01"]),ZvTIL11=pull(get(i)[,"ZvTIL11"]),
                               ZxTIL18=pull(get(i)[,"ZxTIL18"]),ZxTIL25=pull(get(i)[,"ZxTIL25"]))
  }
  if(str_detect(i,"bin2")){
    ALL_collapseddels<-add_row(ALL_collapseddels,
                               CHROM=as.numeric(pull(get(i)[,"CHROM"])),Start=pull(get(i)[,"Start"]),Stop=pull(get(i)[,"Stop"]),
                               REF=pull(get(i)[,"REF"]),ALT=pull(get(i)[,"ALT"]),QUAL=pull(get(i)[,"QUAL"]),M="M2",
                               TdFL=pull(get(i)[,"TdFL"]),TdKS=pull(get(i)[,"TdKS"]),ZdGigi_4to1=pull(get(i)[,"ZdGigi_4to1"]),
                               ZdMomo_4to1=pull(get(i)[,"ZdMomo_4to1"]),ZhRIMHU001=pull(get(i)[,"ZhRIMHU001"]),
                               ZmB73=pull(get(i)[,"ZmB73"]),
                               ZmB97=pull(get(i)[,"ZmB97"]),ZmCML103=pull(get(i)[,"ZmCML103"]),ZmCML228=pull(get(i)[,"ZmCML228"]),
                               ZmCML247=pull(get(i)[,"ZmCML247"]),ZmCML277=pull(get(i)[,"ZmCML277"]),ZmCML322=pull(get(i)[,"ZmCML322"]),
                               ZmCML333=pull(get(i)[,"ZmCML333"]),ZmCML52=pull(get(i)[,"ZmCML52"]),ZmCML69=pull(get(i)[,"ZmCML69"]),
                               ZmHP301=pull(get(i)[,"ZmHP301"]),ZmIL14H=pull(get(i)[,"ZmIL14H"]),ZmKi11=pull(get(i)[,"ZmKi11"]),
                               ZmKi3=pull(get(i)[,"ZmKi3"]),ZmKy21=pull(get(i)[,"ZmKy21"]),ZmM162W=pull(get(i)[,"ZmM162W"]),
                               ZmM37W=pull(get(i)[,"ZmM37W"]), ZmMS71=pull(get(i)[,"ZmMS71"]),ZmMo18W=pull(get(i)[,"ZmMo18W"]),
                               ZmNC350=pull(get(i)[,"ZmNC350"]),ZmNC358=pull(get(i)[,"ZmNC358"]),ZmOh43=pull(get(i)[,"ZmOh43"]),
                               ZmOh7b=pull(get(i)[,"ZmOh7b"]), ZmP39=pull(get(i)[,"ZmP39"]),ZmTx303=pull(get(i)[,"ZmTx303"]),
                               ZmTzi8=pull(get(i)[,"ZmTzi8"]),ZvTIL01=pull(get(i)[,"ZvTIL01"]),ZvTIL11=pull(get(i)[,"ZvTIL11"]),
                               ZxTIL18=pull(get(i)[,"ZxTIL18"]),ZxTIL25=pull(get(i)[,"ZxTIL25"]))
  }
}

ALL_collapseddels<-ALL_collapseddels[-1,] #remove first row of NAs

####Split alt alleles####

#split into different columns
ALL_collapseddels <- ALL_collapseddels %>% mutate(ALT1 = case_when(REF != "collapsedDEL" ~ str_split(ALT, ",",simplify=T)[,1],
                                              REF == "collapsedDEL" ~ "collapsedDEL"),
       ALT2 = case_when(REF != "collapsedDEL" ~ str_split(ALT, ",",simplify=T)[,2],
                        REF == "collapsedDEL" ~ "collapsedDEL"),
       ALT1_length = case_when(str_length(ALT1) > 0  & !ALT1 %in% c("*","collapsedDEL") ~ str_length(ALT1) - str_length(REF),
                               str_length(ALT1) == 0 ~ NA,
                               ALT1 == "collapsedDEL" ~ (Stop - Start)*-1), #*-1 makes it so all Del lengths are < 0 
       ALT2_length = case_when(str_length(ALT2) > 0 & !ALT2 %in% c("*","collapsedDEL") ~ str_length(ALT2) - str_length(REF),
                               str_length(ALT2) == 0 ~ NA,
                               ALT2 == "collapsedDEL" ~ (Stop - Start)*-1), 
       ALT1_type = case_when(ALT1_length < 0 & ALT1 != "collapsedDEL" ~ "del",
                             ALT1_length > 0 & ALT1 != "collapsedDEL"~ "ins",
                             ALT1_length == 0 ~ "snp",
                             is.na(ALT1_length) & !ALT1 %in% c("*","collapsedDEL") ~ NA,
                             ALT1 == "*" ~ "spanning_del",
                             ALT1 == "collapsedDEL" ~ "collapsedDEL"),
       ALT2_type = case_when(ALT2_length < 0 & ALT2 != "collapsedDEL" ~ "del",
                             ALT2_length > 0 & ALT2 != "collapsedDEL" ~ "ins",
                             ALT2_length == 0 ~ "snp",
                             is.na(ALT2_length) & !ALT2 %in% c("*","collapsedDEL") ~ NA,
                             ALT2 == "*" ~ "spanning_del",
                             ALT2 == "collapsedDEL" ~ "collapsedDEL")
)

ALL_collapseddels<-ALL_collapseddels %>% 
  mutate(ALT1_size_class = case_when(ALT1_type %in% c("del","collapsedDEL") & abs(ALT1_length) <= 10 ~ "1-10bp",
                                     ALT1_type %in% c("del","collapsedDEL") & abs(ALT1_length) > 10 & abs(ALT1_length) <= 100~ "11-100bp",
                                     ALT1_type %in% c("del","collapsedDEL") & abs(ALT1_length) > 100 & abs(ALT1_length) <= 1000 ~ "101-1000bp",
                                     ALT1_type %in% c("del","collapsedDEL") & abs(ALT1_length) > 1000 & abs(ALT1_length) <= 10000~ "1001-10Kbp",
                                     ALT1_type %in% c("del","collapsedDEL") & abs(ALT1_length) > 10000 & abs(ALT1_length) <= 100000~ "10K-100Kbp",
                                     ALT1_type %in% c("del","collapsedDEL") & abs(ALT1_length) > 100000 ~ "100K+bp"),
         #only dels (no collapsed dels) for Alt2 because collapsedDels will already be counted in ALT1
         ALT2_size_class = case_when(ALT2_type %in% c("del") & abs(ALT2_length) <= 10 ~ "1-10bp",
                                     ALT2_type %in% c("del") & abs(ALT2_length) > 10 & abs(ALT2_length) <= 100~ "11-100bp",
                                     ALT2_type %in% c("del") & abs(ALT2_length) > 100 & abs(ALT2_length) <= 1000 ~ "101-1000bp",
                                     ALT2_type %in% c("del") & abs(ALT2_length) > 1000 & abs(ALT2_length) <= 10000~ "1001-10Kbp",
                                     ALT2_type %in% c("del") & abs(ALT2_length) > 10000 & abs(ALT2_length) <= 100000~ "10K-100Kbp",
                                     ALT2_type %in% c("del") & abs(ALT2_length) > 100000 ~ "100K+bp")
  ) 
ALL_collapseddels$ALT1_size_class<-ALL_collapseddels$ALT1_size_class %>% factor(levels = c("1-10bp","11-100bp","101-1000bp","1001-10Kbp","10K-100Kbp","100K+bp"))
ALL_collapseddels$ALT2_size_class<-ALL_collapseddels$ALT2_size_class %>% factor(levels = c("1-10bp","11-100bp","101-1000bp","1001-10Kbp","10K-100Kbp","100K+bp"))

ALL_collapseddels %>% group_by(ALT1_size_class, ALT2_size_class) %>% 
  count() %>% na.omit() %>%
  ggplot(aes(x=ALT1_size_class,y=ALT2_size_class, fill = n))+
    geom_tile()+
    theme_bw()+ scale_fill_viridis_c("magma")
ALL_collapseddels %>% group_by(ALT1_type, ALT2_type) %>% 
  count() %>% na.omit() %>%
  ggplot(aes(x=ALT1_type,y=ALT2_type, fill = n))+
  geom_tile()+
  theme_bw()+ scale_fill_viridis_c(option="magma", "Count")+
  xlab("ALT1 Type")+ylab("ALT2 Type")+ 
  scale_x_discrete(labels= c("collapsedDEL"="Collapsed Del.", "spanning_del"="Spanning Del.","del"="Deletion","ins"="Insertion","snp" = "SNP"))+
  scale_y_discrete(labels= c("collapsedDEL"="Collapsed Del.", "spanning_del"="Spanning Del.","del"="Deletion","ins"="Insertion","snp" = "SNP"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/delDescription_type_heatmap.png",
       dpi=300,device="png")

#pivot long so that each row is one allele and its resultant metrics and remove cases where there's no second allele
long_collapsed_dels<-ALL_collapseddels %>% 
  #select(contains("ALT")) %>% 
  pivot_longer(cols = c(ends_with("length")), names_to = "allele_length",names_pattern = "(.*)_length", values_to = "length") %>%
  pivot_longer(cols = c(ends_with("type")), names_to = "allele_type",names_pattern = "(.*)_type", values_to = "type") %>%
  pivot_longer(cols = c(ends_with("size_class")), names_to = "allele_size_class",names_pattern = "(.*)_size_class", values_to = "size_class") %>%
  filter(allele_length == allele_type & allele_length == allele_size_class) %>%
  select(-c(allele_type,allele_size_class)) %>%
  filter(!is.na(type))

colnames(long_collapsed_dels)[45]<-"allele"

ggplot(na.omit(long_collapsed_dels), 
       aes(x=size_class, fill=allele))+
  geom_bar(stat = "count",position = "dodge")+
  facet_wrap(vars(type))+
  theme_bw()+ theme(axis.text.x = element_text(angle = 90))+
  xlab("Size Class")+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/delDescription_sizeclass_barchart.png",
       device="png",dpi=300)

#filter out insertions and SNPs
long_collapsed_dels<-filter(long_collapsed_dels, !(allele %in% "ALT2" & type %in% "collapsedDEL")) 
long_collapsed_dels<-filter(long_collapsed_dels, !(type %in% c("snp","ins"))) 

#These should already be deletions that overlap our reference exons
#However, they are not the deletions that intersect exactly with an exon,
#meaning they may extend beyond the start and stop boundaries of the exon

add_GeneIDs_to_deletions<-function(deletion_df){
  #add columns to the data frame for the exon information
  deletion_df<-add_column(deletion_df, CHROM.ref=NA)
  deletion_df<-add_column(deletion_df,Start.ref=NA)
  deletion_df<-add_column(deletion_df,End.ref=NA)
  deletion_df<-add_column(deletion_df,ID.ref=NA)
  deletion_df<-add_column(deletion_df,Strand.ref=NA)
  deletion_df<-add_column(deletion_df,CDS_ID.ref=NA)
  deletion_df<-add_column(deletion_df,Gene_ID.ref=NA)
  deletion_df<-add_column(deletion_df,CDS_Length.ref=NA)
  for(i in 1:nrow(deletion_df)){
    temp<-filter(ref_Sb313.cds, (as.numeric(CHROM) == deletion_df$CHROM[i] & Start <= deletion_df$Stop[i] & End >= deletion_df$Start[i]))
    deletion_df$CHROM.ref[i]<-paste(temp$CHROM, collapse=":")
    deletion_df$Start.ref[i]<-paste(temp$Start, collapse=":")
    deletion_df$End.ref[i]<-paste(temp$End, collapse=":")
    deletion_df$ID.ref[i]<-paste(temp$ID, collapse=":")
    deletion_df$Strand.ref[i]<-paste(temp$Strand, collapse=":")
    deletion_df$CDS_ID.ref[i]<-paste(temp$CDS_ID, collapse=":")
    deletion_df$Gene_ID.ref[i]<-paste(unique(temp$Gene_ID), collapse=":")
    deletion_df$CDS_Length.ref[i]<-paste(temp$CDS_Length, collapse=":")
  }
  return(deletion_df)
}

ref_Sb313.cds <-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed",
                         col_names = c("CHROM","Start","End","ID","QUAL","Strand"))
ref_Sb313.cds<-ref_Sb313.cds %>% mutate(CDS_ID = str_split(ID, ";",simplify = T)[,1] %>% str_remove_all("ID="),
                                        Gene_ID = str_split(ID, ";",simplify = T)[,2] %>% str_remove_all("Parent="),
                                        CDS_Length = End - Start) %>% select(-QUAL)

#this will take a long time to run
#long_collapsed_dels<-add_GeneIDs_to_deletions(long_collapsed_dels)

#### If adding gene ids to deletions in R (not recommended)####
#add columns to the data frame for the exon information
long_collapsed_dels<-add_column(long_collapsed_dels, CHROM.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,Start.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,End.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,ID.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,Strand.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,CDS_ID.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,Gene_ID.ref=NA)
long_collapsed_dels<-add_column(long_collapsed_dels,CDS_Length.ref=NA)
for(i in 1:nrow(long_collapsed_dels)){
  temp<-filter(ref_Sb313.cds, (as.numeric(CHROM) == long_collapsed_dels$CHROM[i] & Start <= long_collapsed_dels$Stop[i] & End >= long_collapsed_dels$Start[i]))
  long_collapsed_dels$CHROM.ref[i]<-paste(temp$CHROM, collapse=":")
  long_collapsed_dels$Start.ref[i]<-paste(temp$Start, collapse=":")
  long_collapsed_dels$End.ref[i]<-paste(temp$End, collapse=":")
  long_collapsed_dels$ID.ref[i]<-paste(temp$ID, collapse=":")
  long_collapsed_dels$Strand.ref[i]<-paste(temp$Strand, collapse=":")
  long_collapsed_dels$CDS_ID.ref[i]<-paste(temp$CDS_ID, collapse=":")
  long_collapsed_dels$Gene_ID.ref[i]<-paste(unique(temp$Gene_ID), collapse=":")
  long_collapsed_dels$CDS_Length.ref[i]<-paste(temp$CDS_Length, collapse=":")
  if(i %% 100 == 0 ){print(paste((i/nrow(long_collapsed_dels)*100),"% of the way done"))}
  }
#at 30 seconds for 100 rows, that's ~19 hours to finish; 14% in 2 hours suggests it'll finish in 7 hours?
#for troubleshooting:
test<-add_GeneIDs_to_deletions(long_collapsed_dels[1:1000,])
#deletions_with_genes$CHROM.ref<-as.numeric(deletions_with_genes$CHROM.ref)

#This includes all deletions, even those that extend beyond the borders of an exon

####If adding gene ids to deletions in unix (recommended)####

long_collapsed_dels<-mutate(long_collapsed_dels, CHROM=case_when(CHROM == 1 ~ "01",
                                                                 CHROM == 2 ~ "02",
                                                                 CHROM == 3 ~ "03",
                                                                 CHROM == 4 ~ "04",
                                                                 CHROM == 5 ~ "05",
                                                                 CHROM == 6 ~ "06",
                                                                 CHROM == 7 ~ "07",
                                                                 CHROM == 8 ~ "08",
                                                                 CHROM == 9 ~ "09",
                                                                 CHROM == 10 ~ "10"))
write_tsv(long_collapsed_dels, "/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/long_collapsed_dels.tsv")

#in unix
# ml bedtools2
# tail -n +2 /work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/long_collapsed_dels.tsv | bedtools intersect -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed -b - -wo > long_collapsed_dels.withIDs.tsv

long_collapsed_dels<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/long_collapsed_dels.withIDs.tsv",
                              col_names = c("CHROM.ref","Start.ref","Stop.ref","ID","QUAL.ref","Strand.ref",
                                            "CHROM","Start","Stop","REF","ALT",
                                            "QUAL","M","TdFL","TdKS","ZdGigi_4to1",
                                            "ZdMomo_4to1","ZhRIMHU001","ZmB73","ZmB97","ZmCML103",   
                                            "ZmCML228","ZmCML247","ZmCML277","ZmCML322","ZmCML333",   
                                            "ZmCML52","ZmCML69","ZmHP301","ZmIL14H","ZmKi11",    
                                            "ZmKi3","ZmKy21","ZmM162W","ZmM37W","ZmMS71",
                                            "ZmMo18W","ZmNC350","ZmNC358","ZmOh43","ZmOh7b",
                                            "ZmP39","ZmTx303","ZmTzi8","ZvTIL01","ZvTIL11",
                                            "ZxTIL18","ZxTIL25","ALT1","ALT2","allele",
                                            "length","type","size_class","bpOverlap"))

#Where do deletions occur along a gene/cds?

filter(long_collapsed_dels, type == "del") %>% 
  mutate(Start.ref = Start.ref %>% as.numeric(),
         Stop.ref = Stop.ref %>% as.numeric(),
         CDS_Length.ref = Stop.ref - Start.ref,
         position_along_exon=case_when(Strand.ref == "+" ~ (Start - Start.ref)/CDS_Length.ref,
                                       Strand.ref == "-" ~ (Stop.ref - Stop)/CDS_Length.ref)) %>%
  #select(position_along_exon, Start.ref,End.ref, CDS_Length.ref, Start,Stop) %>%
  #filter(position_along_exon >= 0) %>%
  ggplot(aes(x=position_along_exon, y=size_class))+
  geom_jitter(width = 0)+
  geom_violin(alpha=0.5)+
  theme_bw()+ylab("Deletion Size Class")+xlab("Deletion Start Position Along Exon")

#For deletions that are completely within their exons
filter(long_collapsed_dels, type == "del") %>% 
  mutate(Start.ref = Start.ref %>% as.numeric(),
         Stop.ref = Stop.ref %>% as.numeric(),
         CDS_Length.ref = Stop.ref - Start.ref,
         position_along_exon=case_when(Strand.ref == "+" & Start >= Start.ref & Stop <= Stop.ref ~ (Start - Start.ref)/CDS_Length.ref,
                                       Strand.ref == "-" & Start >= Start.ref & Stop <= Stop.ref ~ (Stop.ref - Stop)/CDS_Length.ref)) %>%
  #select(position_along_exon, Start.ref,End.ref, CDS_Length.ref, Start,Stop) %>%
  #filter(position_along_exon >= 0) %>%
  ggplot(aes(x=position_along_exon, y=size_class))+
  geom_jitter(width = 0, aes(color=M))+
  #geom_violin(alpha=0.5, aes(fill = M))+
  theme_bw()+ylab("Deletion Size Class")+xlab("Deletion Start Position Along Exon")

##What we get from this:
#deletion starts or position calls are uniformally distributed across exons, 
#really large deletions tend to start upstream of the exons very far away

#what about at the gene level:
#this will only get the exons with deletions, not the full gene

ref_Sb313.cds<-ref_Sb313.cds %>% group_by(Gene_ID) %>% 
  mutate(Start.gene = min(Start),
         End.gene = max(End),
         Gene_Length.gene = End.gene - Start.gene) %>% ungroup()

long_collapsed_dels %>%   
mutate(Gene_ID.ref = str_split(ID, ";", simplify = T)[,2] %>% str_remove_all("Parent="),
       CDS_ID.ref = str_split(ID, ";", simplify = T)[,1] %>% str_remove_all("ID=")) %>%
  left_join(x=.,y=select(ref_Sb313.cds, CDS_ID, ends_with(".gene")), by=c("CDS_ID.ref"="CDS_ID")) %>%
  mutate(position_along_gene = case_when(Strand.ref == "+" & Start >= Start.gene & Stop <= End.gene ~ (Start - Start.gene)/Gene_Length.gene,
                                         Strand.ref == "-" & Start >= Start.gene & Stop <= End.gene ~ (End.gene - Stop)/Gene_Length.gene)) %>%
  ggplot(aes(x=position_along_gene, y=M))+
  #geom_jitter(width=0, aes(color = size_class))+
  geom_violin(alpha = 0.5)+
  theme_bw()+ylab("Subgenome")+xlab("Deletion Position Along Gene Model")

#This tells us that in M1 there seem to be maybe more deletions that start upstream of the gene? (can't see unless you let the deletions go beyond gene model start and stop)
#and that deletions strictly within the bounds of the gene model occur more at gene model ends than middles
#would likely need a statistical test. 

#Does it change if deletion sizes are kept small
long_collapsed_dels %>% filter(size_class %in% c("1-10bp","11-100bp")) %>%
  mutate(Gene_ID.ref = str_split(ID, ";", simplify = T)[,2] %>% str_remove_all("Parent="),
         CDS_ID.ref = str_split(ID, ";", simplify = T)[,1] %>% str_remove_all("ID=")) %>%
  left_join(x=.,y=select(ungroup(ref_Sb313.cds), CDS_ID, ends_with(".gene")), by=c("CDS_ID.ref"="CDS_ID")) %>%
  mutate(position_along_gene = case_when(Strand.ref == "+" & Start >= Start.gene & Stop <= End.gene ~ (Start - Start.gene)/Gene_Length.gene,
                                         Strand.ref == "-" & Start >= Start.gene & Stop <= End.gene ~ (End.gene - Stop)/Gene_Length.gene)) %>%
  ggplot(aes(x=position_along_gene, y=M))+
  #geom_jitter(width=0, aes(color = size_class))+
  geom_violin(alpha = 0.5)+
  theme_bw()+ylab("Subgenome")+xlab("Deletion Position Along Gene Model")
#No visual changes 

#Which genes/cds don't overlap with a single deletion?

ref_Sb313.cds <-ref_Sb313.cds%>% group_by(Gene_ID) %>% mutate(Total_CDS_Length = sum(CDS_Length))

long_collapsed_dels<-long_collapsed_dels %>%   
  mutate(Gene_ID.ref = str_split(ID, ";", simplify = T)[,2] %>% str_remove_all("Parent="),
         CDS_ID.ref = str_split(ID, ";", simplify = T)[,1] %>% str_remove_all("ID=")) %>%
  left_join(x=.,y=select(ref_Sb313.cds, CDS_ID, ends_with(".gene")), by=c("CDS_ID.ref"="CDS_ID")) %>%
  mutate(position_along_gene = case_when(Strand.ref == "+" & Start >= Start.gene & Stop <= End.gene ~ (Start - Start.gene)/Gene_Length.gene,
                                         Strand.ref == "-" & Start >= Start.gene & Stop <= End.gene ~ (End.gene - Stop)/Gene_Length.gene))

filter(ref_Sb313.cds, !Gene_ID %in% long_collapsed_dels$Gene_ID.ref) %>% write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/NoDelOverlap.ref.genes.tsv")
filter(ref_Sb313.cds, !Gene_ID %in% long_collapsed_dels$Gene_ID.ref) %>% nrow() #6033 rows

filter(ref_Sb313.cds, !CDS_ID %in% long_collapsed_dels$CDS_ID.ref)%>% write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/NoDelOverlap.ref.cds.tsv")

filter(ref_Sb313.cds, Gene_ID %in% long_collapsed_dels$Gene_ID.ref) %>% write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/DelOverlap.ref.genes.tsv")
filter(ref_Sb313.cds, Gene_ID %in% long_collapsed_dels$Gene_ID.ref) %>% nrow() #63294 rows

filter(ref_Sb313.cds, CDS_ID %in% long_collapsed_dels$CDS_ID.ref)%>% write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/DelOverlap.ref.cds.tsv")


#### Calculate the number of bp in a deletion ####
long_collapsed_dels %>% group_by(Gene_ID, M) %>% summarize(bpDelOverlap.gene = sum(bpOverlap,na.rm=T)) %>%
  left_join(y=ref_Sb313.cds, x=., by="Gene_ID") %>% 
  mutate(PercentDelOverlap.gene = (bpDelOverlap.gene/Gene_Length.gene)*100) %>%
  #ggplot(aes(x=Gene_Length.gene, y=PercentDelOverlap.gene))+
  ggplot(aes(x=Gene_Length.gene, y=bpDelOverlap.gene))+
  geom_point(aes(color = M))+
  geom_abline(slope = 1)+
  theme_bw()
#What this shows is that the small genes have lots over different deletions or are encountering extremely large deletions
#Or there are polymorphic deletions, such that for a given genome, only part of the total number of bp that overlap with a deletion are actually deleted

#Need to calculate it based by genome

CDS_bpDel_byGenomeM<-long_collapsed_dels %>%  filter(type != "spanning_del")%>%
  mutate(TdFL.delbp = case_when(TdFL==0 ~ 0, 
                                TdFL==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                TdFL==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                TdFL > 0 & type == "collapsedDEL"~ bpOverlap),
         TdKS.delbp = case_when(TdKS==0 ~ 0, 
                                TdKS==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                TdKS==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                TdKS > 0 & type == "collapsedDEL"~ bpOverlap),
         ZdGigi_4to1.delbp = case_when(ZdGigi_4to1==0 ~ 0, 
                                       ZdGigi_4to1==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                       ZdGigi_4to1==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                       ZdGigi_4to1 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZdMomo_4to1.delbp = case_when(ZdMomo_4to1==0 ~ 0, 
                                       ZdMomo_4to1==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                       ZdMomo_4to1==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                       ZdMomo_4to1 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZhRIMHU001.delbp = case_when(ZhRIMHU001==0 ~ 0, 
                                      ZhRIMHU001==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                      ZhRIMHU001==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                      ZhRIMHU001 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmB73.delbp = case_when(ZmB73==0 ~ 0, 
                                 ZmB73==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmB73==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmB73 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmB97.delbp = case_when(ZmB97==0 ~ 0, 
                                 ZmB97==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmB97==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmB97 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML103.delbp = case_when(ZmCML103==0 ~ 0, 
                                    ZmCML103==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML103==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML103 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML228.delbp = case_when(ZmCML228==0 ~ 0, 
                                    ZmCML228==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML228==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML228 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML247.delbp = case_when(ZmCML247==0 ~ 0, 
                                    ZmCML247==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML247==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML247 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML277.delbp = case_when(ZmCML277==0 ~ 0, 
                                    ZmCML277==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML277==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML277 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML322.delbp = case_when(ZmCML322==0 ~ 0, 
                                    ZmCML322==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML322==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML322 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML333.delbp = case_when(ZmCML333==0 ~ 0, 
                                    ZmCML333==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML333==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML333 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML52.delbp = case_when(ZmCML52==0 ~ 0, 
                                   ZmCML52==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmCML52==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmCML52 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML69.delbp = case_when(ZmCML69==0 ~ 0, 
                                   ZmCML69==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmCML69==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmCML69 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmHP301.delbp = case_when(ZmHP301==0 ~ 0, 
                                   ZmHP301==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmHP301==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmHP301 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmIL14H.delbp = case_when(ZmIL14H==0 ~ 0, 
                                   ZmIL14H==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmIL14H==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmIL14H > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmKi11.delbp = case_when(ZmKi11==0 ~ 0, 
                                  ZmKi11==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmKi11==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmKi11 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmKi3.delbp = case_when(ZmKi3==0 ~ 0, 
                                 ZmKi3==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmKi3==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmKi3 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmKy21.delbp = case_when(ZmKy21==0 ~ 0, 
                                  ZmKy21==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmKy21==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmKy21 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmM162W.delbp = case_when(ZmM162W==0 ~ 0, 
                                   ZmM162W==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmM162W==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmM162W > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmM37W.delbp = case_when(ZmM37W==0 ~ 0, 
                                  ZmM37W==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmM37W==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmM37W > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmMS71.delbp = case_when(ZmMS71==0 ~ 0, 
                                  ZmMS71==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmMS71==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmMS71 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmMo18W.delbp = case_when(ZmMo18W==0 ~ 0, 
                                   ZmMo18W==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmMo18W==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmMo18W > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmNC350.delbp = case_when(ZmNC350==0 ~ 0, 
                                   ZmNC350==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmNC350==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmNC350 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmNC358.delbp = case_when(ZmNC358==0 ~ 0, 
                                   ZmNC358==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmNC358==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmNC358 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmOh43.delbp = case_when(ZmOh43==0 ~ 0, 
                                  ZmOh43==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmOh43==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmOh43 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmOh7b.delbp = case_when(ZmOh7b==0 ~ 0, 
                                  ZmOh7b==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmOh7b==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmOh7b > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmP39.delbp = case_when(ZmP39==0 ~ 0, 
                                 ZmP39==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmP39==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmP39 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmTx303.delbp = case_when(ZmTx303==0 ~ 0, 
                                   ZmTx303==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmTx303==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmTx303 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmTzi8.delbp = case_when(ZmTzi8==0 ~ 0, 
                                  ZmTzi8==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmTzi8==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmTzi8 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZvTIL01.delbp = case_when(ZvTIL01==0 ~ 0, 
                                   ZvTIL01==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZvTIL01==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZvTIL01 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZvTIL11.delbp = case_when(ZvTIL11==0 ~ 0, 
                                   ZvTIL11==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZvTIL11==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZvTIL11 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZxTIL18.delbp = case_when(ZxTIL18==0 ~ 0, 
                                   ZxTIL18==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZxTIL18==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZxTIL18 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZxTIL25.delbp = case_when(ZxTIL25==0 ~ 0, 
                                   ZxTIL25==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZxTIL25==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZxTIL25 > 0 & type == "collapsedDEL"~ bpOverlap),) %>% 
    select(CDS_ID.ref,Gene_ID.ref,M,bpOverlap, ends_with(".delbp"))%>%
  group_by(CDS_ID.ref, M) %>% summarize(across(ends_with(".delbp"), \(x) sum(x,na.rm =T)))

#Are there any cases where the bpOverlap of a deletion is greater than the length of the gene? 
filter(long_collapsed_dels, bpOverlap > Gene_Length.gene) %>% nrow()

CDS_bpDel_byGenomeM<-left_join(y=ref_Sb313.cds, x=CDS_bpDel_byGenomeM, by=c("CDS_ID.ref"="CDS_ID"))  
  #mutate(across(ends_with(".delbp"), ~ (.x/CDS_Length))*100) #this is for converting to proportion

CDS_bpDel_byGenomeM %>% 
  pivot_longer(ends_with(".delbp"), names_to = "Genome", values_to = "Delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp")) %>%
  ggplot(aes(x=Genome, y=Delbp))+
  geom_violin(aes(fill = M))+
  coord_flip()+
  theme_bw()+ ylab("Deleted bps in an exon")

CDS_bpDel_byGenomeM %>% 
  mutate(across(ends_with(".delbp"), ~ (.x/CDS_Length))) %>%
  pivot_longer(ends_with(".delbp"), names_to = "Genome", values_to = "Delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp")) %>%
  ggplot(aes(x=Genome, y=Delbp))+
  geom_violin(aes(fill = M))+
  coord_flip()+
  theme_bw()+ ylab("Proportion of Exon in Deleted bp")

#Plot cds length by bp del but for each genome
CDS_bpDel_byGenomeM %>% 
  pivot_longer(cols = ends_with(".delbp"),names_to = "Genome", values_to = "Delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp")) %>%
  ggplot(aes(x=CDS_Length, y=Delbp, color = Genome))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1)+
  theme_bw()+
  xlab("Exon bp length")+ylab("Deleted bp of Exon")+
  scale_color_manual(values = genome_colors)+
  facet_wrap(vars(M))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/CDSdeletedBPsByCDSLength.Genome.png", device="png",dpi = 300, width=6, height = 5)

gene_bpDel_byGenomeM<-long_collapsed_dels %>%  filter(type != "spanning_del")%>%
  mutate(TdFL.delbp = case_when(TdFL==0 ~ 0, 
                                TdFL==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                TdFL==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                TdFL > 0 & type == "collapsedDEL"~ bpOverlap),
         TdKS.delbp = case_when(TdKS==0 ~ 0, 
                                TdKS==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                TdKS==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                TdKS > 0 & type == "collapsedDEL"~ bpOverlap),
         ZdGigi_4to1.delbp = case_when(ZdGigi_4to1==0 ~ 0, 
                                       ZdGigi_4to1==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                       ZdGigi_4to1==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                       ZdGigi_4to1 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZdMomo_4to1.delbp = case_when(ZdMomo_4to1==0 ~ 0, 
                                       ZdMomo_4to1==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                       ZdMomo_4to1==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                       ZdMomo_4to1 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZhRIMHU001.delbp = case_when(ZhRIMHU001==0 ~ 0, 
                                      ZhRIMHU001==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                      ZhRIMHU001==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                      ZhRIMHU001 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmB73.delbp = case_when(ZmB73==0 ~ 0, 
                                 ZmB73==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmB73==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmB73 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmB97.delbp = case_when(ZmB97==0 ~ 0, 
                                 ZmB97==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmB97==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmB97 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML103.delbp = case_when(ZmCML103==0 ~ 0, 
                                    ZmCML103==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML103==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML103 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML228.delbp = case_when(ZmCML228==0 ~ 0, 
                                    ZmCML228==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML228==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML228 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML247.delbp = case_when(ZmCML247==0 ~ 0, 
                                    ZmCML247==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML247==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML247 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML277.delbp = case_when(ZmCML277==0 ~ 0, 
                                    ZmCML277==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML277==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML277 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML322.delbp = case_when(ZmCML322==0 ~ 0, 
                                    ZmCML322==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML322==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML322 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML333.delbp = case_when(ZmCML333==0 ~ 0, 
                                    ZmCML333==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                    ZmCML333==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                    ZmCML333 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML52.delbp = case_when(ZmCML52==0 ~ 0, 
                                   ZmCML52==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmCML52==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmCML52 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmCML69.delbp = case_when(ZmCML69==0 ~ 0, 
                                   ZmCML69==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmCML69==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmCML69 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmHP301.delbp = case_when(ZmHP301==0 ~ 0, 
                                   ZmHP301==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmHP301==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmHP301 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmIL14H.delbp = case_when(ZmIL14H==0 ~ 0, 
                                   ZmIL14H==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmIL14H==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmIL14H > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmKi11.delbp = case_when(ZmKi11==0 ~ 0, 
                                  ZmKi11==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmKi11==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmKi11 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmKi3.delbp = case_when(ZmKi3==0 ~ 0, 
                                 ZmKi3==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmKi3==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmKi3 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmKy21.delbp = case_when(ZmKy21==0 ~ 0, 
                                  ZmKy21==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmKy21==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmKy21 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmM162W.delbp = case_when(ZmM162W==0 ~ 0, 
                                   ZmM162W==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmM162W==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmM162W > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmM37W.delbp = case_when(ZmM37W==0 ~ 0, 
                                  ZmM37W==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmM37W==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmM37W > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmMS71.delbp = case_when(ZmMS71==0 ~ 0, 
                                  ZmMS71==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmMS71==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmMS71 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmMo18W.delbp = case_when(ZmMo18W==0 ~ 0, 
                                   ZmMo18W==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmMo18W==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmMo18W > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmNC350.delbp = case_when(ZmNC350==0 ~ 0, 
                                   ZmNC350==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmNC350==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmNC350 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmNC358.delbp = case_when(ZmNC358==0 ~ 0, 
                                   ZmNC358==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmNC358==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmNC358 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmOh43.delbp = case_when(ZmOh43==0 ~ 0, 
                                  ZmOh43==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmOh43==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmOh43 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmOh7b.delbp = case_when(ZmOh7b==0 ~ 0, 
                                  ZmOh7b==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmOh7b==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmOh7b > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmP39.delbp = case_when(ZmP39==0 ~ 0, 
                                 ZmP39==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                 ZmP39==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                 ZmP39 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmTx303.delbp = case_when(ZmTx303==0 ~ 0, 
                                   ZmTx303==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZmTx303==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZmTx303 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZmTzi8.delbp = case_when(ZmTzi8==0 ~ 0, 
                                  ZmTzi8==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                  ZmTzi8==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                  ZmTzi8 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZvTIL01.delbp = case_when(ZvTIL01==0 ~ 0, 
                                   ZvTIL01==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZvTIL01==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZvTIL01 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZvTIL11.delbp = case_when(ZvTIL11==0 ~ 0, 
                                   ZvTIL11==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZvTIL11==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZvTIL11 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZxTIL18.delbp = case_when(ZxTIL18==0 ~ 0, 
                                   ZxTIL18==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZxTIL18==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZxTIL18 > 0 & type == "collapsedDEL"~ bpOverlap),
         ZxTIL25.delbp = case_when(ZxTIL25==0 ~ 0, 
                                   ZxTIL25==1 & allele == "ALT1" & type == "del" ~ bpOverlap, 
                                   ZxTIL25==2 & allele =="ALT2" & type == "del"~ bpOverlap,
                                   ZxTIL25 > 0 & type == "collapsedDEL"~ bpOverlap),) %>% 
  select(Gene_ID.ref,M,bpOverlap, ends_with(".delbp"))%>%
  group_by(Gene_ID.ref, M) %>% summarize(across(ends_with(".delbp"), \(x) sum(x,na.rm =T))) %>%
  unique()

gene_bpDel_byGenomeM<-left_join(y=unique(select(ref_Sb313.cds, Gene_ID, Start.gene, End.gene, Gene_Length.gene)), 
                                x=gene_bpDel_byGenomeM, by=c("Gene_ID.ref"="Gene_ID"))  

gene_bpDel_byGenomeM %>% 
  pivot_longer(ends_with(".delbp"), names_to = "Genome", values_to = "Delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp")) %>%
  ggplot(aes(x=Genome, y=Delbp))+
  geom_violin(aes(fill = M))+
  coord_flip()+
  theme_bw()+ ylab("Deleted bps in a gene")

gene_bpDel_byGenomeM %>% 
  mutate(across(ends_with(".delbp"), ~ (.x/Gene_Length.gene))) %>%
  pivot_longer(ends_with(".delbp"), names_to = "Genome", values_to = "Delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp")) %>%
  ggplot(aes(x=Genome, y=Delbp))+
  geom_violin(aes(fill = M))+
  coord_flip()+
  theme_bw()+ ylab("Proportion of Gene in Deleted bp")

#Plot cds length by bp del but for each genome
gene_bpDel_byGenomeM %>% 
  pivot_longer(cols = ends_with(".delbp"),names_to = "Genome", values_to = "Delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp")) %>%
  ggplot(aes(x=Gene_Length.gene, y=Delbp, color = Genome))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1)+
  theme_bw()+
  xlab("Gene bp length")+ylab("Deleted bp of Gene")+
  scale_color_manual(values = genome_colors)+
  facet_wrap(vars(M))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/genedeletedBPsByGeneLength.Genome.png", device="png",dpi = 300, width=6, height = 5)

####Compare fractionation status from the GVCFs to VCF calls####
#if a gene or cds id is  the long_collapsed_dels list, it would be called fractionated
#but a given genome might not have that deletion... So use the CDS_bpDel_byGenomeM
colnames(CDS_bpDel_byGenomeM)
colnames(long.full.fractionation.status)

CDS_bpDel_byGenomeM<-pivot_longer(CDS_bpDel_byGenomeM, ends_with("delbp"), names_to = "Genome", values_to = "delbp") %>%
  mutate(Genome = str_remove_all(Genome, ".delbp"))

gvcf_vcf_agreement<-mutate(long.full.fractionation.status, CDS_ID = str_split(ID, ";",simplify = T)[,1] %>% str_remove_all("ID=")) %>%
  inner_join(x=., y=CDS_bpDel_byGenomeM, by=c("Genome"="Genome","CDS_ID"="CDS_ID.ref")) %>%
  mutate(agreement = case_when(M == "M1" & delbp > 0 & Status %in% c("M2_Retained","Both_Lost","M1_Lost:M2_NA") ~ "Agree",
                               M == "M1" & delbp == 0 & Status %in% c("Both_Retained","M1_Retained","M1_Retained:M2_NA") ~ "Agree",
                               M == "M1" & Status %in% c("M1_NA:M2_Lost","M1_NA:M2_Retained","Both_NA") ~ "Unclear",
                               M == "M1" & delbp > 0 & Status %in% c("Both_Retained","M1_Retained","M1_Retained:M2_NA") ~ "Disagree",
                               M == "M1" & delbp == 0 & Status %in% c("M2_Retained","Both_Lost","M1_Lost:M2_NA") ~ "Disagree",
                               
                               M == "M2" & delbp > 0 & Status %in% c("M1_Retained","Both_Lost","M1_NA:M2_Lost") ~ "Agree",
                               M == "M2" & delbp == 0 & Status %in% c("Both_Retained","M2_Retained","M1_NA:M2_Retained") ~ "Agree",
                               M == "M2" & Status %in% c("M1_Lost:M2_NA","M1_Retained:M2_NA","Both_NA") ~ "Unclear",
                               M == "M2" & delbp > 0 & Status %in% c("Both_Retained","M2_Retained","M1_NA:M2_Retained") ~ "Disagree",
                               M == "M2" & delbp == 0 & Status %in% c("M1_Retained","Both_Lost","M1_NA:M2_Lost")  ~ "Disagree")
  )

gvcf_vcf_agreement %>% group_by(Status, M, agreement) %>% count()
filter(gvcf_vcf_agreement, agreement == "Disagree") %>% View()

gvcf_vcf_agreement %>% group_by(Status, M, Genome, agreement) %>% count() %>%
  ggplot(aes(x=Genome, y=Status))+
  geom_tile(aes(fill=n))+
  facet_grid(rows = vars(agreement),cols = vars(M), scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

#let's try a wide pivot because I think the difference in status calling is tripping things up
#but filling values with 0 might also be an issue for non-alignment

#because we're only looking at deletions in VCF (not retention)
#agreement can only be about fractionation call (both lost, 1 or other lost)
CDS_bpDel_byGenomeM<- CDS_bpDel_byGenomeM %>% mutate(value = T) %>%
  pivot_wider(names_from = M) %>%
  pivot_longer(cols=ends_with(".delbp"), names_to = "Genome", values_to = "Delbp")  %>%
  mutate(Genome = str_remove_all(Genome,".delbp"))

CDS_bpDel_byGenomeM<-mutate(CDS_bpDel_byGenomeM,
         Status = case_when(Delbp > 0 & M1 == T & is.na(M2) ~ "M1_Lost",
                            Delbp > 0 & M2 == T & is.na(M1) ~ "M2_Lost",
                            Delbp == 0 & M1 == T & is.na(M2) ~ "M1_Ret",
                            Delbp == 0 & M2 == T & is.na(M1) ~ "M2_Ret",
                            Delbp > 0 & M1 == T & M2 == T ~ "Both_Lost",
                            Delbp == 0 & M2 == T & M1 == T ~ "Both_Ret",
                            .default = "Unknown"))

gvcf_vcf_agreement<-mutate(long.full.fractionation.status, CDS_ID = str_split(ID, ";",simplify = T)[,1] %>% str_remove_all("ID=")) %>%
  inner_join(x=., y=CDS_bpDel_byGenomeM, by=c("Genome"="Genome","CDS_ID"="CDS_ID.ref"))

gvcf_vcf_agreement <-gvcf_vcf_agreement %>% 
  mutate(Agreement = case_when(
    Status.x == "M1_Lost:M2_NA" & Status.y == "M1_Lost" ~ "Agree",
    Status.x == "M1_Retained:M2_NA" & Status.y==  "M1_Ret" ~ "Agree",
    Status.x == "Both_Retained" & Status.y== "Both_Ret" ~ "Agree",
    Status.x == "M2_Retained" & Status.y== "M2_Ret" ~ "Agree",
    Status.x == "M1_Retained" & Status.y== "M1_Ret" ~ "Agree",
    Status.x == "M1_NA:M2_Lost" & Status.y== "M2_Lost" ~ "Agree",
    Status.x == "M1_NA:M2_Retained" & Status.y == "M2_Lost" ~ "Agree",
    Status.x == "Both_NA" ~ "Uncertain",
    Status.x == "M1_Lost:M2_NA" & Status.y  %in% c("M2_Lost", "M2_Ret", "Both_Lost","Both_Ret") ~ "Uncertain",
    Status.x == "M1_Retained:M2_NA" & Status.y %in% c("M2_Lost", "M2_Ret", "Both_Lost","Both_Ret") ~ "Uncertain",
    Status.x == "M1_NA:M2_Lost" & Status.y %in% c("M1_Lost", "M1_Ret", "Both_Lost","Both_Ret") ~ "Uncertain",
    Status.x == "M1_NA:M2_Retained" & Status.y %in% c("M2_Lost", "M2_Ret", "Both_Lost","Both_Ret") ~ "Uncertain",
    .default = "Disagree"
  ))

#to look at counts
gvcf_vcf_agreement %>% group_by(Agreement) %>% count()
gvcf_vcf_agreement %>% filter(Agreement == "Disagree") %>% group_by(Status.x, Status.y) %>% count() %>% arrange(-n) %>%print(n=25)

gvcf_vcf_agreement %>%  group_by(Status.x, Status.y) %>% count() %>% 
  mutate(Percentage = (n/nrow(gvcf_vcf_agreement))*100) %>%
  ggplot(aes(x=Status.x, y=Status.y))+
  geom_tile(aes(fill = Percentage)) + 
  theme_bw()+xlab("GVCF Calls")+ylab("VCF Calls")+ theme(axis.text.x = element_text(angle = 90))+
  scale_fill_viridis_c()+
  scale_x_discrete(labels=c("M1_Lost:M2_NA" = "M1 Lost, M2 NA","M1_Retained:M2_NA" ="M1 Ret., M2 NA",
                            "Both_Lost"="Both Lost","Both_Retained"="Both Ret.","M2_Retained"="M2 Ret., M1 Lost",
                            'M1_Retained'="M1 Ret., M2 Lost","M1_NA:M2_Lost"="M1 NA, M2 Lost","M1_NA:M2_Retained"="M1 NA, M2 Ret.", "Both_NA"="Both NA"))+
  scale_y_discrete(labels = c("M1_Lost" = "M1 Lost","M2_Lost"="M2 Lost","M1_Ret"="M1 Ret.","M2_Ret"="M2 Ret.","Both_Ret"="Both Ret.","Both_Lost"="Both Lost"))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/gvcf_vs_vcf_agreement.heatmap.png",
       device="png",dpi = 300)

#Overall they mostly agree though it's hard to compare because fractionation by current VCF analysis is only keeping deletions, and thereby doesn't report retention

#### Does multiple deletions help explain paraphyly?####
M1_paraphyly<-filter(full.fractionation.status, M1.Timing.NoZnZd == "paraphyly") %>% select(contains("ID"), ends_with(".M1"))
M1_paraphyly <- M1_paraphyly %>% select(-c(contains("ZnPI"),"ZdGigi.M1","ZdMomo.M1"))
M1_paraphyly_dels<-long_collapsed_dels %>% filter(CDS_ID.ref %in% M1_paraphyly$CDS_ID & M == "M1")

M1_paraphyly[1,] %>% pivot_longer(cols = ends_with("M1")) %>% select(name,value) %>% print(n=35)
filter(M1_paraphyly_dels, CDS_ID.ref == "Sobic.001G002000.1.v3.1.CDS.7") %>% View()

filter(gvcf_vcf_agreement, CDS_ID == "Sobic.001G002000.1.v3.1.CDS.7") %>% select(Genome, starts_with("Stat"),Agreement)%>% print(n=35)

#Perhaps a better question is: is there more disagreement between methods for paraphyly calls?

p_call<-M1_paraphyly[1,] %>% pivot_longer(cols = ends_with("M1")) %>% mutate(name = str_remove_all(name, ".M1"))
a_call<-filter(gvcf_vcf_agreement, CDS_ID == p_call$CDS_ID[1])
#how to find the genomes causing the paraphyly?
#Make the assumption that genomes that have the opposite call to the majority of genomes are causing the paraphyly
p_list<-filter(p_call, value != names(which.max(table(p_call$value)))) %>% select(name)
#get the number of disagreements for genomes in the paraphyly causing list
filter(a_call, Genome %in% p_list$name & Agreement == "Disagree" & M1 == T) %>% nrow()
#vs the total number of disagreements for that CDS total
filter(a_call, Agreement == "Disagree") %>% nrow()
#then see if those genomes == "disagree" in a_call

paraphyly_agreement<-tibble(CDS_ID=NA, Total_Disagreements=NA, SusParaphylyDisagreements=NA, Total_SusParaphyly=NA)
for(i in 1:nrow(M1_paraphyly)){
  p_call<-M1_paraphyly[i,] %>% pivot_longer(cols = ends_with("M1")) %>% mutate(name = str_remove_all(name, ".M1"))
  a_call<-filter(gvcf_vcf_agreement, CDS_ID == p_call$CDS_ID[1] & M1 == T)
  paraphyly_agreement<-add_row(paraphyly_agreement,
                               CDS_ID = p_call$CDS_ID[1],
                               Total_Disagreements=filter(a_call, Agreement == "Disagree") %>% nrow(), 
                               SusParaphylyDisagreements=filter(a_call, Genome %in% pull(select(filter(p_call, value != names(which.max(table(p_call$value)))), name)) & Agreement == "Disagree") %>% nrow(),
                               Total_SusParaphyly=filter(p_call, value != names(which.max(table(p_call$value)))) %>% nrow()
  )
}
paraphyly_agreement<-paraphyly_agreement[-1,]

paraphyly_agreement %>% 
  ggplot(aes(x=Total_SusParaphyly))+
  geom_histogram(binwidth = 1)
#what we would expect given that previously we'd seen a couple genomes deviate from the rest for fractionation call being the most common paraphyly patterns
paraphyly_agreement %>%
  mutate(Per_totalDisagree = (Total_Disagreements/35)*100,
         Per_SusParaphylyDisagreements = (SusParaphylyDisagreements/Total_SusParaphyly)*100) %>%
  ggplot()+
  geom_histogram(aes(x=Per_totalDisagree),binwidth = 5)+
  geom_histogram(aes(x=Per_SusParaphylyDisagreements),binwidth = 5,fill = "darkblue",alpha =0.5)
#What this tells me is that CDS paraphly is likely related to how things are called
#given that a lot of it is in the 100% disagrees, but that's a hard metric to trust
#because vcf vs. gvcf calls aren't a truly equal comparison to make and so it's hard 
#to differentiate true signal from this noise
#Trying to look at this with genes... that'd be the convergence analysis




mutate(PercentDelOverlap.gene = (bpDelOverlap.gene/Gene_Length.gene)*100) %>%
  #ggplot(aes(x=Gene_Length.gene, y=PercentDelOverlap.gene))+
  ggplot(aes(x=Gene_Length.gene, y=bpDelOverlap.gene))+
  geom_point(aes(color = M))+
  geom_abline(slope = 1)+
  theme_bw()




#just working out the basics using TdFL
temp<-filter(long_collapsed_dels, str_detect(Gene_ID.ref, ref_Sb313.cds$Gene_ID[322])) #good test indices are 12, 322
temp<-mutate(temp, 
             TdFL.del.bp = case_when(TdFL == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                     TdFL == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                     TdFL == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
             TdKS.del.bp = case_when(TdKS == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                     TdKS == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                     TdKS == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
             ZmCML228.del.bp = case_when(ZmCML228 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmCML228 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmCML228 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length))
#when multiple genes are involved
#if Start < min(Start.ref) of a gene, set to start.ref 
#if Stop > max(End.ref) of a gene, set to End.ref 
#this includes intron bp deleted similar to the more simple cases
#length is stop - start after reset
for(i in 1:nrow(temp)){ #go through each row of temp
  if(str_detect(temp[i,"Gene_ID.ref"],":")){ #if the row has multiple gene ID
    #figure out which sections belong to that gene
    indices<-str_split(temp$CDS_ID.ref[i], ":",simplify = T) %>% str_detect(ref_Sb313.cds$Gene_ID[322]) %>% which()
    start<-str_split(temp$Start.ref[i],":",simplify=T)[,indices] %>% as.numeric() %>% min() %>% c(temp$Start[i]) %>% max()
    end.bp<-str_split(temp$End.ref[i],":",simplify=T)[,indices] %>% as.numeric() %>% max() %>% c(temp$Stop[i]) %>% min()
    #so long as it hasn't been assigned del.bp and the genotype is not NA

    if(temp$type[i] == "collapsedDEL"){ #if the type is collapsedDEL
      if(is.na(temp$TdFL.del.bp[i]) & temp$TdFL[i] > 0){temp$TdFL.del.bp[i]<-end.bp-start}
      if(is.na(temp$TdKS.del.bp[i]) & temp$TdKS[i] > 0){temp$TdKS.del.bp[i]<-end.bp-start}
      if(is.na(temp$ZmCML228.del.bp[i]) & temp$ZmCML228[i] > 0){temp$ZmCML228.del.bp[i]<-end.bp-start}
    }else{ #if not collapsed DEL nor a spanning deletion
      if(temp$type[i] != "spanning_del"){ #then match allele and genotype before assigning end-start
        if((temp$allele[i] == "ALT1" & temp$TdFL[i] == 1) | (temp$allele[i] == "ALT2" & temp$TdFL[i] == 2)){temp$TdFL.del.bp[i]<-end.bp-start}
        if((temp$allele[i] == "ALT1" & temp$TdKS[i] == 1) | (temp$allele[i] == "ALT2" & temp$TdKS[i] == 2)){temp$TdKSL.del.bp[i]<-end.bp-start}
        if((temp$allele[i] == "ALT1" & temp$ZmCML228[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML228[i] == 2)){temp$ZmCML228.del.bp[i]<-end.bp-start}
      }
    }
  }
}

temp %>% select(TdFL, TdKS, ZmCML228, type,allele, length, ends_with(".del.bp"))

temp %>% group_by(M) %>% summarize(across(str_subset(colnames(temp), ".del.bp"),\(x) sum(abs(x),na.rm = T)))                                     


#make it into a function:
calculate_total_del_bp<-function(delDF, refGeneID){
  temp<-filter(delDF, str_detect(Gene_ID.ref, refGeneID)) #good test indices are 12, 322
  if(nrow(temp) > 0){
  temp<-mutate(temp, 
               TdFL.del.bp = case_when(TdFL == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                       TdFL == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                       TdFL == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               TdKS.del.bp = case_when(TdKS == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                       TdKS == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                       TdKS == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZdGigi_4to1.del.bp = case_when(ZdGigi_4to1 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                              ZdGigi_4to1 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                              ZdGigi_4to1 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZdMomo_4to1.del.bp = case_when(ZdMomo_4to1 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                              ZdMomo_4to1 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                              ZdMomo_4to1 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZhRIMHU001.del.bp = case_when(ZhRIMHU001 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                             ZhRIMHU001 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                             ZhRIMHU001 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZvTIL01.del.bp = case_when(ZvTIL01 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZvTIL01 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZvTIL01 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZvTIL11.del.bp = case_when(ZvTIL11 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZvTIL11 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZvTIL11 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZxTIL18.del.bp = case_when(ZxTIL18 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZxTIL18 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZxTIL18 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZxTIL25.del.bp = case_when(ZxTIL25 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZxTIL25 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZxTIL25 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmB73.del.bp = case_when(ZmB73 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                        ZmB73 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                        ZmB73 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmB97.del.bp = case_when(ZmB97 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                        ZmB97 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                        ZmB97 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML103.del.bp = case_when(ZmCML103 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                           ZmCML103 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                           ZmCML103 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML228.del.bp = case_when(ZmCML228 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                           ZmCML228 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                           ZmCML228 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML247.del.bp = case_when(ZmCML247 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                           ZmCML247 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                           ZmCML247 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML277.del.bp = case_when(ZmCML277 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                           ZmCML277 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                           ZmCML277 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML322.del.bp = case_when(ZmCML322 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                           ZmCML322 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                           ZmCML322 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML333.del.bp = case_when(ZmCML333 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                           ZmCML333 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                           ZmCML333 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML52.del.bp = case_when(ZmCML52 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmCML52 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmCML52 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmCML69.del.bp = case_when(ZmCML69 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmCML69 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmCML69 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmHP301.del.bp = case_when(ZmHP301 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmHP301 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmHP301 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmIL14H.del.bp = case_when(ZmIL14H == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmIL14H == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmIL14H == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmKi11.del.bp = case_when(ZmKi11 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmKi11 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmKi11 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmKi3.del.bp = case_when(ZmKi3 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                        ZmKi3 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                        ZmKi3 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmKy21.del.bp = case_when(ZmKy21 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmKy21 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmKy21 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmM162W.del.bp = case_when(ZmM162W == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmM162W == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmM162W == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmM37W.del.bp = case_when(ZmM37W == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmM37W == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmM37W == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmMo18W.del.bp = case_when(ZmMo18W == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmMo18W == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmMo18W == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmMS71.del.bp = case_when(ZmMS71 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmMS71 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmMS71 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmNC350.del.bp = case_when(ZmNC350 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmNC350 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmNC350 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmNC358.del.bp = case_when(ZmNC358 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmNC358 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmNC358 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmOh43.del.bp = case_when(ZmOh43 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmOh43 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmOh43 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmOh7b.del.bp = case_when(ZmOh7b == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmOh7b == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmOh7b == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmP39.del.bp = case_when(ZmP39 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                        ZmP39 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                        ZmP39 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmTx303.del.bp = case_when(ZmTx303 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                          ZmTx303 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                          ZmTx303 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length),
               ZmTzi8.del.bp = case_when(ZmTzi8 == 0 ~ 0, #simple cases, spanning deletions will be NA (keeps bp from being counted twice)
                                         ZmTzi8 == 1 & allele == "ALT1" & str_detect(Gene_ID.ref,":",negate = T) ~ length,
                                         ZmTzi8 == 2 & allele == "ALT2" & str_detect(Gene_ID.ref,":",negate = T) ~ length))
  #when multiple genes are involved
  #if Start < min(Start.ref) of a gene, set to start.ref 
  #if Stop > max(End.ref) of a gene, set to End.ref 
  #this includes intron bp deleted similar to the more simple cases
  #length is stop - start after reset
  for(i in 1:nrow(temp)){ #go through each row of temp
    if(str_detect(temp[i,"Gene_ID.ref"],":")){ #if the row has multiple gene ID
      #figure out which sections belong to that gene
      indices<-str_split(temp$CDS_ID.ref[i], ":",simplify = T) %>% str_detect(refGeneID) %>% which()
      start<-str_split(temp$Start.ref[i],":",simplify=T)[,indices] %>% as.numeric() %>% min() %>% c(temp$Start[i]) %>% max()
      end.bp<-str_split(temp$End.ref[i],":",simplify=T)[,indices] %>% as.numeric() %>% max() %>% c(temp$Stop[i]) %>% min()
      #so long as it hasn't been assigned del.bp and the genotype is not NA
      if(temp$type[i] == "collapsedDEL"){ #if the type is collapsedDEL
        if(is.na(temp$TdFL.del.bp[i]) & temp$TdFL[i] > 0){temp$TdFL.del.bp[i]<-end.bp-start}
        if(is.na(temp$TdKS.del.bp[i]) & temp$TdKS[i] > 0){temp$TdKS.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZdGigi_4to1.del.bp[i]) & temp$ZdGigi_4to1[i] > 0){temp$ZdGigi_4to1.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZdMomo_4to1.del.bp[i]) & temp$ZdMomo_4to1[i] > 0){temp$ZdMomo_4to1.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZhRIMHU001.del.bp[i]) & temp$ZhRIMHU001[i] > 0){temp$ZhRIMHU001.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZvTIL01.del.bp[i]) & temp$ZvTIL01[i] > 0){temp$ZvTIL01.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZvTIL11.del.bp[i]) & temp$ZvTIL11[i] > 0){temp$ZvTIL11.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZxTIL18.del.bp[i]) & temp$ZxTIL18[i] > 0){temp$ZxTIL18.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZxTIL25.del.bp[i]) & temp$ZxTIL25[i] > 0){temp$ZxTIL25.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmB73.del.bp[i]) & temp$ZmB73[i] > 0){temp$ZmB73.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmB97.del.bp[i]) & temp$ZmB97[i] > 0){temp$ZmB97.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML103.del.bp[i]) & temp$ZmCML103[i] > 0){temp$ZmCML103.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML228.del.bp[i]) & temp$ZmCML228[i] > 0){temp$ZmCML228.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML247.del.bp[i]) & temp$ZmCML247[i] > 0){temp$ZmCML247.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML277.del.bp[i]) & temp$ZmCML277[i] > 0){temp$ZmCML277.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML322.del.bp[i]) & temp$ZmCML322[i] > 0){temp$ZmCML322.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML333.del.bp[i]) & temp$ZmCML333[i] > 0){temp$ZmCML333.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML52.del.bp[i]) & temp$ZmCML52[i] > 0){temp$ZmCML52.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmCML69.del.bp[i]) & temp$ZmCML69[i] > 0){temp$ZmCML69.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmHP301.del.bp[i]) & temp$ZmHP301[i] > 0){temp$ZmHP301.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmIL14H.del.bp[i]) & temp$ZmIL14H[i] > 0){temp$ZmIL14H.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmKi11.del.bp[i]) & temp$ZmKi11[i] > 0){temp$ZmKi11.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmKi3.del.bp[i]) & temp$ZmKi3[i] > 0){temp$ZmKi3.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmKy21.del.bp[i]) & temp$ZmKy21[i] > 0){temp$ZmKy21.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmM162W.del.bp[i]) & temp$ZmM162W[i] > 0){temp$ZmM162W.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmM37W.del.bp[i]) & temp$ZmM37W[i] > 0){temp$ZmM37W.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmMo18W.del.bp[i]) & temp$ZmMo18W[i] > 0){temp$ZmMo18W.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmMS71.del.bp[i]) & temp$ZmMS71[i] > 0){temp$ZmMS71.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmNC350.del.bp[i]) & temp$ZmNC350[i] > 0){temp$ZmNC350.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmNC358.del.bp[i]) & temp$ZmNC358[i] > 0){temp$ZmNC358.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmOh43.del.bp[i]) & temp$ZmOh43[i] > 0){temp$ZmOh43.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmOh7b.del.bp[i]) & temp$ZmOh7b[i] > 0){temp$ZmOh7b.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmP39.del.bp[i]) & temp$ZmP39[i] > 0){temp$ZmP39.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmTx303.del.bp[i]) & temp$ZmTx303[i] > 0){temp$ZmTx303.del.bp[i]<-end.bp-start}
        if(is.na(temp$ZmTzi8.del.bp[i]) & temp$ZmTzi8[i] > 0){temp$ZmTzi8.del.bp[i]<-end.bp-start}
      }else{#if not collapsed DEL nor a spanning deletion
        if(temp$type[i] != "spanning_del"){ #then match allele and genotype before assigning end-start
          if((temp$allele[i] == "ALT1" & temp$TdFL[i] == 1) | (temp$allele[i] == "ALT2" & temp$TdFL[i] == 2)){temp$TdFL.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$TdKS[i] == 1) | (temp$allele[i] == "ALT2" & temp$TdKS[i] == 2)){temp$TdKS.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZdGigi_4to1[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZdGigi_4to1[i] == 2)){temp$ZdGigi_4to1.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZdMomo_4to1[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZdMomo_4to1[i] == 2)){temp$ZdMomo_4to1.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZhRIMHU001[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZhRIMHU001[i] == 2)){temp$ZhRIMHU001.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZvTIL01[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZvTIL01[i] == 2)){temp$ZvTIL01.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZvTIL11[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZvTIL11[i] == 2)){temp$ZvTIL11.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZxTIL18[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZxTIL18[i] == 2)){temp$ZxTIL18.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZxTIL25[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZxTIL25[i] == 2)){temp$ZxTIL25.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmB73[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmB73[i] == 2)){temp$ZmB73.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmB97[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmB97[i] == 2)){temp$ZmB97.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML103[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML103[i] == 2)){temp$ZmCML103.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML228[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML228[i] == 2)){temp$ZmCML228.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML247[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML247[i] == 2)){temp$ZmCML247.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML277[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML277[i] == 2)){temp$ZmCML277.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML322[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML322[i] == 2)){temp$ZmCML322.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML333[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML333[i] == 2)){temp$ZmCML333.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML52[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML52[i] == 2)){temp$ZmCML52.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmCML69[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmCML69[i] == 2)){temp$ZmCML69.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmHP301[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmHP301[i] == 2)){temp$ZmHP301.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmIL14H[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmIL14H[i] == 2)){temp$ZmIL14H.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmKi11[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmKi11[i] == 2)){temp$ZmKi11.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmKi3[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmKi3[i] == 2)){temp$ZmKi3.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmKy21[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmKy21[i] == 2)){temp$ZmKy21.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmM162W[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmM162W[i] == 2)){temp$ZmM162W.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmM37W[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmM37W[i] == 2)){temp$ZmM37W.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmMo18W[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmMo18W[i] == 2)){temp$ZmMo18W.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmMS71[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmMS71[i] == 2)){temp$ZmMS71.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmNC350[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmNC350[i] == 2)){temp$ZmNC350.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmNC358[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmNC358[i] == 2)){temp$ZmNC358.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmOh43[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmOh43[i] == 2)){temp$ZmOh43.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmOh7b[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmOh7b[i] == 2)){temp$ZmOh7b.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmP39[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmP39[i] == 2)){temp$ZmP39.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmTx303[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmTx303[i] == 2)){temp$ZmTx303.del.bp[i]<-end.bp-start}
          if((temp$allele[i] == "ALT1" & temp$ZmTzi8[i] == 1) | (temp$allele[i] == "ALT2" & temp$ZmTzi8[i] == 2)){temp$ZmTzi8.del.bp[i]<-end.bp-start}
        }
      }
    }
  }
  }
  return(temp %>% group_by(M) %>% summarize(across(str_subset(colnames(temp), ".del.bp"),\(x) sum(abs(x),na.rm = T))))
}

calculate_total_del_bp(test, ref_Sb313.cds$Gene_ID[322]) %>% mutate(across(!M, \(x) (x/ref_Sb313.cds$Length.gene[322])*100)) %>% view()
#note that you can still end up with more than 100% of the bp lost, if deletions overlap each other but weren't collapsed
calculate_total_del_bp(test, ref_Sb313.cds$Gene_ID[322]) %>% mutate_all(\(x) (x/ref_Sb313.cds$Length.gene[322])*100) %>% str()

del_bp_by_gene<-tibble(Gene_ID=NA, M= NA, TdFL=NA,TdKS=NA,ZdGigi_4to1=NA,ZdMomo_4to1=NA,ZhRIMHU001=NA,ZvTIL01=NA,ZvTIL11=NA,
                           ZxTIL18=NA,ZxTIL25=NA,ZmB73=NA,ZmB97=NA,ZmCML103=NA,ZmCML228=NA,ZmCML247=NA,ZmCML277=NA,ZmCML322=NA,
                           ZmCML333=NA,ZmCML52=NA,ZmCML69=NA,ZmHP301=NA,ZmIL14H=NA,ZmKi11=NA,ZmKi3=NA,ZmKy21=NA,ZmM162W=NA,
                           ZmM37W=NA,ZmMo18W=NA,ZmMS71=NA,ZmNC350=NA,ZmNC358=NA,ZmOh43=NA,ZmOh7b=NA,ZmP39=NA,ZmTx303=NA,ZmTzi8=NA)
for(i in unique(ref_Sb313.cds$Gene_ID)){
  temp<-calculate_total_del_bp(test, i) #%>% mutate_all(\(x) (x/ref_Sb313.cds$Length.gene[i])*100)
  if(!all(is.na(temp[1,]))){del_bp_by_gene<-add_row(del_bp_by_gene, Gene_ID = i,M=temp$M,
                                                        TdFL=temp$TdFL.del.bp,TdKS=temp$TdKS.del.bp,ZdGigi_4to1=temp$ZdGigi_4to1.del.bp,
                                                        ZdMomo_4to1=temp$ZdMomo_4to1.del.bp,ZhRIMHU001=temp$ZhRIMHU001.del.bp,ZvTIL01=temp$ZvTIL01.del.bp,
                                                        ZvTIL11=temp$ZvTIL11.del.bp,
                                                        ZxTIL18=temp$ZxTIL18.del.bp,ZxTIL25=temp$ZxTIL25.del.bp,ZmB73=temp$ZmB73.del.bp,
                                                        ZmB97=temp$ZmB97.del.bp,ZmCML103=temp$ZmCML103.del.bp,ZmCML228=temp$ZmCML228.del.bp,
                                                        ZmCML247=temp$ZmCML247.del.bp,ZmCML277=temp$ZmCML277.del.bp,ZmCML322=temp$ZmCML322.del.bp,
                                                        ZmCML333=temp$ZmCML333.del.bp,ZmCML52=temp$ZmCML52.del.bp,ZmCML69=temp$ZmCML69.del.bp,
                                                        ZmHP301=temp$ZmHP301.del.bp,ZmIL14H=temp$ZmIL14H.del.bp,ZmKi11=temp$ZmKi11.del.bp,
                                                        ZmKi3=temp$ZmKi3.del.bp,ZmKy21=temp$ZmKy21.del.bp,ZmM162W=temp$ZmM162W.del.bp,
                                                        ZmM37W=temp$ZmM37W.del.bp,ZmMo18W=temp$ZmMo18W.del.bp,ZmMS71=temp$ZmMS71.del.bp,
                                                        ZmNC350=temp$ZmNC350.del.bp,
                                                        ZmNC358=temp$ZmNC358.del.bp,ZmOh43=temp$ZmOh43.del.bp,ZmOh7b=temp$ZmOh7b.del.bp,
                                                        ZmP39=temp$ZmP39.del.bp,ZmTx303=temp$ZmTx303.del.bp,ZmTzi8=temp$ZmTzi8.del.bp)}
}
del_bp_by_gene<-del_bp_by_gene[-1,] #removes first row of all NAs

del_bp_by_gene<-select(ref_Sb313.cds, Gene_ID, Gene_Length.gene) %>% unique() %>% inner_join(x=.,y=del_bp_by_gene,by="Gene_ID") %>%
  mutate(across(!c(matches("Gene"),M), \(x) (x/Gene_Length.gene)*100))

pivot_longer(del_bp_by_gene, cols = !c(matches("Gene"),M), names_to = "Genome",values_to = "Percent_BP_del") %>%
  filter(Percent_BP_del <= 100) %>%
  ggplot(aes(x=Percent_BP_del,y=M,fill=Genome, color=Genome))+
  geom_jitter(width = 0)#geom_violin()

pivot_longer(del_bp_by_gene, cols = !c(matches("Gene"),M), names_to = "Genome",values_to = "Percent_BP_del") %>%
  group_by(Genome) %>% 
  summarize(min.bp=min(Percent_BP_del),
            max.bp=max(Percent_BP_del),
            mean.bp=mean(Percent_BP_del),
            median.bp=median(Percent_BP_del)) %>% view()
#but this doesn't solve the issue if there's a large deletion that spans part of two genes
#because then both genes would have the total deletion length applied to them...







#since the objects should still be by chromosome
#filter out intergenic only indels
##so long as the end of the deletion is after the start of the exon
##and the start of the deletion is before the end of the exon
##then the deletion should intersect with the exon

findGenicDels<-function(obj, refchr){
  temp_ref<-filter(ref_Sb313.cds, CHROM == refchr)
  newdf<-tibble(CHROM = NA, POS=NA,End=NA, REF=NA, ALT=NA,QUAL=NA,ALT1_Length=NA,ALT2_Length=NA,
                ALT1_type=NA, ALT2_type=NA,
                TdFL=NA,TdKS=NA,ZdGigi_4to1=NA,ZdMomo_4to1=NA,
                ZhRIMHU001=NA,ZmB73=NA,ZmB97=NA,ZmCML103=NA,ZmCML228=NA,ZmCML247=NA,ZmCML277=NA,
                ZmCML322=NA,ZmCML333=NA,ZmCML52=NA,ZmCML69=NA,ZmHP301=NA,ZmIL14H=NA,ZmKi11=NA,
                ZmKi3=NA,ZmKy21=NA,ZmM162W=NA,ZmM37W=NA,ZmMS71=NA,ZmMo18W=NA,
                ZmNC358=NA,ZmNC350=NA,ZmOh43=NA,ZmOh7b=NA,ZmP39=NA,ZmTx303=NA,ZmTzi8=NA,
                ZvTIL01=NA,ZvTIL11=NA,ZxTIL18=NA,ZxTIL25=NA,ZnPI615697_4to1=NA,
                RefExon.Start=NA, RefExon.End=NA, RefExon.ID=NA, RefExon.Strand=NA,RefExon.CDS_ID=NA,
                RefExon.Gene_ID=NA, RefExon.CDS_Length=NA)
  for(i in 1:nrow(temp_ref)){
    temp_query<-filter(obj, CHROM == refchr & POS <= temp_ref$End[i] & END >= temp_ref$Start[i])
    if(nrow(temp_query) > 0){
      newdf<-add_row(newdf,CHROM = pull(temp_query[,1]), POS=pull(temp_query[,2]),End=pull(temp_query[,"END"]), 
                     REF=pull(temp_query[,4]), ALT=pull(temp_query[,5]),QUAL=pull(temp_query[,6]),ALT1_Length=pull(temp_query[,48]),
                     ALT2_Length=pull(temp_query[,49]),ALT1_type=pull(temp_query[,50]),ALT2_type=pull(temp_query[,51]),
                     TdFL=pull(temp_query[,10]),TdKS=pull(temp_query[,11]), 
                     ZdGigi_4to1=pull(temp_query[,12]),ZdMomo_4to1=pull(temp_query[,13]),
                     ZhRIMHU001=pull(temp_query[,14]),ZmB73=pull(temp_query[,15]),ZmB97=pull(temp_query[,16]),
                     ZmCML103=pull(temp_query[,17]),
                     ZmCML228=pull(temp_query[,18]),ZmCML247=pull(temp_query[,19]),ZmCML277=pull(temp_query[,20]),
                     ZmCML322=pull(temp_query[,21]),ZmCML333=pull(temp_query[,22]),ZmCML52=pull(temp_query[,23]),ZmCML69=pull(temp_query[,24]),
                     ZmHP301=pull(temp_query[,25]),ZmIL14H=pull(temp_query[,26]),ZmKi11=pull(temp_query[,27]),
                     ZmKi3=pull(temp_query[,28]),
                     ZmKy21=pull(temp_query[,29]),ZmM162W=pull(temp_query[,30]),ZmM37W=pull(temp_query[,31]),
                     ZmMS71=pull(temp_query[,32]),ZmMo18W=pull(temp_query[,33]),
                     ZmNC358=pull(temp_query[,35]),ZmNC350=pull(temp_query[,34]),ZmOh43=pull(temp_query[,36]),
                     ZmOh7b=pull(temp_query[,37]),ZmP39=pull(temp_query[,38]),ZmTx303=pull(temp_query[,39]),ZmTzi8=pull(temp_query[,40]),
                     ZvTIL01=pull(temp_query[,42]),ZvTIL11=pull(temp_query[,43]),ZxTIL18=pull(temp_query[,44]),ZxTIL25=pull(temp_query[,45]),
                     ZnPI615697_4to1=pull(temp_query[,41]),
                     RefExon.Start=temp_ref$Start[i], RefExon.End=temp_ref$End[i], RefExon.ID=temp_ref$ID[i], 
                     RefExon.Strand=temp_ref$Strand[i],RefExon.CDS_ID=temp_ref$CDS_ID[i],
                     RefExon.Gene_ID=temp_ref$Gene_ID[i], RefExon.CDS_Length=temp_ref$CDS_Length[i])
    }
  }
  return(newdf)
}

#test.genic.M1<-head(test.M1, n=5000) %>% findGenicDels(., "02")
test.genic.M1<-findGenicDels(test.M1, "02")
test.genic.M1<-test.genic.M1[-1,] #removes initializing row
#changes the periods in genotypes to "NA"
test.genic.M1<-mutate(test.genic.M1, across(c(starts_with("Z"),starts_with("T")),~ str_replace(string = .x,pattern = "\\.",replacement = "NA")))
#changes the "NA" to NA
test.genic.M1[test.genic.M1=="NA"]<-NA

####Deletions biallelic (only deletion within exon?) vs multi-allelic (multiple deletions within an exon?)####
###defined by variant with the same starting POS for biallelic
delbyexon.counts<-test.genic.M1 %>% group_by(RefExon.CDS_ID) %>% count()

test.biallelic.M1<-filter(test.genic.M1, RefExon.CDS_ID %in% pull(select(filter(delbyexon.counts, n == 1), "RefExon.CDS_ID")))
test.multiallelic.M1<-filter(test.genic.M1, RefExon.CDS_ID %in% pull(select(filter(delbyexon.counts, n > 1), "RefExon.CDS_ID")))

#nrow of these dfs would give the number of deletions under each
#uniq(RefExon.CDS_ID) %>% length() would give the number of unique exons under each

####How many are nested (deletions physically overlap)?####
##This would be spanning deletions
#number of spanning deletions
filter(test.genic.M1, ALT1_type == "spanning_del" | ALT2_type == "spanning_del") %>% nrow() 
#number of unique exons that include a spanning deletion variant
filter(test.genic.M1, ALT1_type == "spanning_del" | ALT2_type == "spanning_del") %>% select(RefExon.CDS_ID) %>% unique() %>% nrow()

####How many deletions extend past exon boundries vs. confined within?####

#deletions within exon boundaries
#will include spanning deletions if the other alt is not a spanning deletion but is confined within boundaries
test.genic.M1 %>% filter(POS > RefExon.Start & End < RefExon.End)

####How many confined deletions are in-frame vs frame shifting?####
#in-frame
test.genic.M1 %>% filter(POS > RefExon.Start & End < RefExon.End) %>% 
  filter(ALT1_Length %% 3 == 0 | ALT2_Length %% 3 == 0)
#frame-shifting
test.genic.M1 %>% filter(POS > RefExon.Start & End < RefExon.End) %>% 
  filter(ALT1_Length %% 3 != 0 | ALT2_Length %% 3 != 0)

########################################################
####Redoing timing estimates using the combined vcf####

#For each variant row: 
##split into those that have REF, ALT1, ALT2 (if applicable) and skip those with NAs
##For those genomes that share ALT1 or ALT2
###is it a monophyletic group?
####if yes: assign to a node
####if no: paraphyly #This would be suggestive of ILS, whereas multiple polymorphic would be multiple origins?
######So this would be assigning deletions individually

#Splitting and redoing genotypes such that ALT1 and ALT2 are on separate lines and spanning del genotypes = NA
test.genic.long.M1<-pivot_longer(test.genic.M1, cols = c(starts_with("Z"),starts_with("T")),
             names_to = "Genome",
             values_to = "Genotype") 
test.genic.long.M1 <-test.genic.long.M1%>% mutate(Genotype = case_when(Genotype == 1 & ALT1_type == "spanning_del" ~ NA,
                                                   Genotype == 2 & ALT2_type == "spanning_del" ~ NA,
                                                   .default = Genotype))

test.genic.long.M1<-test.genic.long.M1 %>% pivot_wider(names_from = Genome, values_from = Genotype) 


test.genic.long.M1<-test.genic.long.M1 %>% pivot_longer(cols = c(contains("type")), 
                               names_to = "ALT_type_name",
                               values_to = "ALT_type",
                               values_drop_na = T) %>% 
  pivot_longer(cols = c(ALT1_Length, ALT2_Length),
               names_to = "ALT_length_name",
               values_to = "ALT_length",
               values_drop_na = T) 
test.genic.long.M1<-filter(test.genic.long.M1,ALT_type != "spanning_del")
test.genic.long.M1<-filter(test.genic.long.M1, (ALT_type_name == "ALT1_type" & ALT_length_name == "ALT1_Length") | (ALT_type_name == "ALT2_type" & ALT_length_name == "ALT2_Length"))
test.genic.long.M1<-filter(test.genic.long.M1, ALT_type == "del")

####Create a revised timing categories that can work even with 1's and 2's####
timing_categories.revised<-tibble(Node = c("basal_lost","Tripsacum_lost","Zea_lost",
                                           "Diplo_lost","Mays_lost",#"huehuetenangensis_lost",
                                           "MexParvMaize_lost", "Mex_lost","ParvMaize_lost",
                                           "Parv_lost","Maize_lost","TropMaize_lost","TempMaize_lost",
                                           "Trop_woTzi8_lost","Trop_NCs_lost","Trop_CMLsKis_lost",
                                           "Trop_CMLsKisWoCML333_lost","Trop_CMLSKisWoCML333orCML322",
                                           "Trop_CMLsalone_lost","Trop_CMLsandKisRecent_lost", "Trop_CML247and277_lost",
                                           "Trop_CML69and52_lost","Trop_KisandCML228_lost","Trop_Ki3and228_lost",
                                           "Temp_woM162W_lost","Temp_OhandKy21","Temp_FlintBsandMS71","Temp_Ohonly",
                                           "Temp_flintonly","Temp_BsandMS71","Temp_sweetonly","temp_B97andMS71",
                                           tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]),
                                  Genomes = c(paste(tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)],collapse = ":"),#all genomes
                                              paste(tripsacinae_genome_IDs[c(1:2)],collapse = ":"),#just trip
                                              paste(tripsacinae_genome_IDs[c(3:4,8:24,27:32,34,36:39)],collapse = ":"),#just zea
                                              paste(tripsacinae_genome_IDs[c(3:4)],collapse = ":"),#just diplo
                                              paste(tripsacinae_genome_IDs[c(8:24,27:32,34,36:39)],collapse = ":"),#just mays
                                              #paste(tripsacinae_genome_IDs[c(8)],collapse = ":"), #just huehuetenagensis
                                              paste(tripsacinae_genome_IDs[c(9:24,27:32,34,36:39)],collapse = ":"),#mex, parv, maize
                                              paste(tripsacinae_genome_IDs[c(38:39)],collapse = ":"),#just mex
                                              paste(tripsacinae_genome_IDs[c(9:24,27:32,34,36:37)],collapse = ":"), #parv, maize
                                              paste(tripsacinae_genome_IDs[c(36:37)],collapse = ":"),#just parv
                                              paste(tripsacinae_genome_IDs[c(9:24,27:32,34)],collapse = ":"),#just maize
                                              paste(tripsacinae_genome_IDs[c(11:18,21:22,28:29,34)],collapse = ":"),#just tropical
                                              paste(tripsacinae_genome_IDs[c(9:10,19:20,23:24,27,30:32)],collapse = ":"),#just temperate
                                              paste(tripsacinae_genome_IDs[c(11:18,21:22,28:29)],collapse = ":"),#Trop-no Tzi8
                                              paste(tripsacinae_genome_IDs[c(28:29)],collapse = ":"),#Trop-just NCs
                                              paste(tripsacinae_genome_IDs[c(11:18,21:22)],collapse = ":"),#Trop-CMLs and Kis
                                              paste(tripsacinae_genome_IDs[c(11:15,17:18,21:22)],collapse = ":"),#Trop CMLs and Kis without CML333
                                              paste(tripsacinae_genome_IDs[c(11:14,17:18,21:22)],collapse = ":"),#Trop CMLs and Kis without CML322 and CML333
                                              paste(tripsacinae_genome_IDs[c(11,13:14)],collapse = ":"),#Trop CMLs alone (CML103,247,277)
                                              paste(tripsacinae_genome_IDs[c(12,17:18,21:22)],collapse = ":"),#Trop CMLs and Kis (CML52,69,228, and Ki3 and 11)
                                              paste(tripsacinae_genome_IDs[c(13:14)],collapse = ":"), #Trop CML247 and CML277
                                              paste(tripsacinae_genome_IDs[c(17:18)],collapse = ":"), #Trop CML69 and 52
                                              paste(tripsacinae_genome_IDs[c(12,21:22)],collapse = ":"),#Trop Ki3 CML228 and Ki11
                                              paste(tripsacinae_genome_IDs[c(12,22)],collapse = ":"),#TropKi3 and CML228
                                              paste(tripsacinae_genome_IDs[c(9:10,19:20,23,27,30:32)],collapse = ":"),#temperate without M162W
                                              paste(tripsacinae_genome_IDs[c(23,30:31)],collapse = ":"),#Temparate Oh and Ky21
                                              paste(tripsacinae_genome_IDs[c(9:10,19:20,27,32)],collapse = ":"),#temparate Flints, Bs and MS71
                                              paste(tripsacinae_genome_IDs[c(30:31)],collapse = ":"),#temperate Oh only
                                              paste(tripsacinae_genome_IDs[c(19:20,32)],collapse = ":"),#temparate Flints only
                                              paste(tripsacinae_genome_IDs[c(9:10,27)],collapse = ":"),#temparate B73 B97 MS71
                                              paste(tripsacinae_genome_IDs[c(20,32)],collapse = ":"),#temparate sweet only
                                              paste(tripsacinae_genome_IDs[c(10,27)],collapse = ":"),#temperate B97 and MS71
                                              tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)] #each individual genome
                                              )
                                  )
#add columns for each genome
for(i in tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]){timing_categories.revised[,i]<-NA}
#add T/F for the categories by genomes
for(i in tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]){
  for(j in 1:nrow(timing_categories.revised)){
    if(str_detect(string = timing_categories.revised$Genomes[j], pattern = i) == T){
      timing_categories.revised[j,i]<-T
    }else{
      timing_categories.revised[j,i]<-F
    }
  }
}

##Adapt previous timing estimate loop to the above dataframe:

for(i in 1:nrow(test.genic.long.M1)){
  if(all(!is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]]))){ #if all have a value, no na
    if(test.genic.long.M1$ALT_type_name[i] == "ALT1_type"){ #if row is ALT1
      #creates the string of genomes in the same order and format as the strings above
      pat<-which(colSums(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]] == 1, na.rm = TRUE) > 0) %>% names() %>% paste(collapse = ":")
      node<-filter(timing_categories.revised, Genomes == pat)
      if(nrow(node) > 0){
        test.genic.long.M1[i,"M1.TimingNode.NoZnZd"]<-node$Node[1]
      }else{test.genic.long.M1[i,"M1.TimingNode.NoZnZd"]<-"paraphyly"}
    }
    if(test.genic.long.M1$ALT_type_name[i] == "ALT2_type"){ #if row is ALT2
      #creates the string of genomes in the same order and format as the strings above
      pat<-which(colSums(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]] == 2, na.rm = TRUE) > 0) %>% names() %>% paste(collapse = ":")
      node<-filter(timing_categories.revised, Genomes == pat)
      if(nrow(node) > 0){
        test.genic.long.M1[i,"M1.TimingNode.NoZnZd"]<-node$Node[1]
      }else{test.genic.long.M1[i,"M1.TimingNode.NoZnZd"]<-"paraphyly"}
    }
  }
  if(any(is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]]))){
    test.genic.long.M1[i,"M1.TimingNode.NoZnZd"]<-"Unaligned"
  }
}
###counts
test.genic.long.M1 %>% group_by(M1.TimingNode.NoZnZd) %>% count()
test.genic.long.M1 %>% group_by(M1.TimingNode.NoZnZd) %>% count() %>% 
  filter(M1.TimingNode.NoZnZd != "Unaligned") %>% ungroup()%>%select(n) %>% colSums()
#6758 total deletions that aren't in unaligned

###plot timing by deletion
test.genic.long.M1 %>% filter(M1.TimingNode.NoZnZd != "Unaligned") %>%
  group_by(M1.TimingNode.NoZnZd) %>% count() %>%
  mutate(Percent = (n/6758)*100) %>%
ggplot(aes(x=M1.TimingNode.NoZnZd, y=Percent))+
  geom_bar(stat="identity")+ coord_flip()+
  ggtitle("Timing with no NA genotypes allowed")

####What are the levels of unalignment among deletions####
aligned_del_count<-0
unaligned_del<-c() #number of genomes per deletion that are unaligned
for(i in 1:nrow(test.genic.long.M1)){
  if(all(!is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]]))){ #if all have a value, no na
  aligned_del_count<-aligned_del_count+1
  }else{#if at least one of those columns has an na
    unaligned_del<-c(unaligned_del, rowSums(is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]])))
  }
}
as_tibble(unaligned_del) %>%
  ggplot(aes(x=value))+geom_histogram(bins = 32)

#count of all unaligned deletions (total)
sum(unaligned_del)#128822
as_tibble(unaligned_del) %>% group_by(value) %>% count() %>% mutate(Percentage = (n/128822)*100) %>%
  ggplot(aes(x=value, y=Percentage))+geom_bar(stat="identity")

####Try timing assignment while allowing for up to 2 NAs
test.genic.long.M1$TimingNode.w2NAs<-NA
for(i in 1:nrow(test.genic.long.M1)){
  if(all(!is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]])) | rowSums(is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]])) <= 2){ #if all have a value or up to 2 na's
    if(test.genic.long.M1$ALT_type_name[i] == "ALT1_type"){ #if row is ALT1
      #creates the string of genomes in the same order and format as the strings above
      pat<-which(colSums(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]] == 1, na.rm = TRUE) > 0) %>% names() %>% paste(collapse = ":")
      node<-filter(timing_categories.revised, Genomes == pat)
      if(nrow(node) > 0){
        test.genic.long.M1[i,"TimingNode.w2NAs"]<-node$Node[1]
      }else{test.genic.long.M1[i,"TimingNode.w2NAs"]<-"paraphyly"}
    }
    if(test.genic.long.M1$ALT_type_name[i] == "ALT2_type"){ #if row is ALT2
      #creates the string of genomes in the same order and format as the strings above
      pat<-which(colSums(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]] == 2, na.rm = TRUE) > 0) %>% names() %>% paste(collapse = ":")
      node<-filter(timing_categories.revised, Genomes == pat)
      if(nrow(node) > 0){
        test.genic.long.M1[i,"TimingNode.w2NAs"]<-node$Node[1]
      }else{test.genic.long.M1[i,"TimingNode.w2NAs"]<-"paraphyly"}
    }
  }
  if(rowSums(is.na(test.genic.long.M1[i,tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]])) > 2){
    test.genic.long.M1[i,"TimingNode.w2NAs"]<-"Unaligned"
  }
}

test.genic.long.M1 %>% group_by(TimingNode.w2NAs) %>% count() %>% 
  filter(TimingNode.w2NAs != "Unaligned") %>% ungroup()%>%select(n) %>% colSums()
#14375 total deletions that aren't in unaligned

###plot timing by deletion
test.genic.long.M1 %>% filter(TimingNode.w2NAs != "Unaligned") %>%
  group_by(TimingNode.w2NAs) %>% count() %>%
  mutate(Percent = (n/14375)*100) %>%
  ggplot(aes(x=TimingNode.w2NAs, y=Percent))+
  geom_bar(stat="identity")+ coord_flip() +
  ggtitle("Timing allowing 2 genomes to be NA")
