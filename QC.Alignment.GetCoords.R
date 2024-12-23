library(tidyverse)

#Read in anchors
for(i in tripsacinae_genome_IDs){
  assign(paste0(i,".anchors"), read.table(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/Sb313_",i,"_anchorwave.anchors"), header = TRUE))
}

#In bash make a gene bed file that could be imported here
#grep "gene" Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 | awk -v OFS='\t' '{print $1,$4,$5,$9,$7}' - > Sb313.gene.bed
Sb313.gene<-read_tsv(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.gene.bed", col_names = c("CHROM","Start","End","ID","Strand"))
Sb313.gene<-mutate(Sb313.gene, 
                   gene_length = End - Start,
                   Gene_ID = str_extract_all(ID, "ID.+?v3.1", simplify = T) %>% str_remove_all("ID=") %>% str_remove_all(".v3.1"))

#Add in the reference coordinates
#Going to start off with all the genes used in the analysis
#in case a large number aren't in anchors
QC.genes<-filter(Sb313.gene, str_detect(Gene_ID, paste(str_remove(ref_Sb313.cds$Gene_ID, pattern = ".[0-9].v3.[0-9]"), collapse = "|"))) 

#Need to be able to separate the query coordinates by subgenome (because some Sb genes will have 2 sets of query coordinates)
for(i in tripsacinae_genome_IDs){
  assign(paste0("parsingCoords_",i), read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/parsing_coordinates/",i,".bins.txt"), col_names = c("queryChr","RefChr","Bin")))
}
subgenome_key<-tibble(Genome = NA, Subgenome = NA, queryChr=NA, RefChr=NA)
for(i in tripsacinae_genome_IDs){
  subgenome_key<-add_row(subgenome_key, Genome = i, Subgenome = pull(get(paste0("parsingCoords_",i))[,3]),
                         queryChr = pull(get(paste0("parsingCoords_",i))[,1]),
                         RefChr = pull(get(paste0("parsingCoords_",i))[,2]))
}
subgenome_key<-subgenome_key[-1,]
subgenome_key$Subgenome<-subgenome_key$Subgenome %>% str_replace_all(pattern="bin",replacement = "M")

#make it 2x 
QC.genes<-QC.genes %>% add_row(CHROM = QC.genes$CHROM, Start = QC.genes$Start, End=QC.genes$End, ID=QC.genes$ID, Strand=QC.genes$Strand, gene_length = QC.genes$gene_length, Gene_ID=QC.genes$Gene_ID)

QC.genes <-add_column(QC.genes, Subgenome = c(rep("M1", 0.5*nrow(QC.genes)), rep("M2",0.5*nrow(QC.genes)))) #create a M1 and M2 copy of each Sorghum gene

#Add in Query coordinates (practice with TdFL)
TdFL.anchors %>% head()
df<-filter(TdFL.anchors, str_detect(gene, paste(QC.genes$Gene_ID, collapse = "|"))) %>% select(c(contains("query"),"strand","gene")) %>% mutate(gene = gene %>% str_remove_all(".[0-9].v3.1")) %>%
  left_join(filter(subgenome_key, Genome == "TdFL"), by="queryChr")

for(t in "TdFL"){
  QC.genes[,paste0(t,"_chr")]<-NA
  QC.genes[,paste0(t,"_start")]<-NA
  QC.genes[,paste0(t,"_end")]<-NA
  QC.genes[,paste0(t,"_strand")]<-NA
}

for(i in 1:nrow(QC.genes)){
  parsing<-filter(subgenome_key, Genome == "TdFL" & RefChr == QC.genes$CHROM[i] & Subgenome == QC.genes$Subgenome[i])
  if(nrow(filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i]& RefChr == QC.genes$CHROM[i])) >0){
    QC.genes[i,paste0(t,"_chr")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,1]
    QC.genes[i,paste0(t,"_start")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,2]
    QC.genes[i,paste0(t,"_end")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,3]
    QC.genes[i,paste0(t,"_strand")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,4]
  }
}

#Make it a loop
for(t in 1:length(tripsacinae_genome_IDs)){
  df<-filter(get(paste0(tripsacinae_genome_IDs[t],".anchors")),str_detect(gene, paste(QC.genes$Gene_ID, collapse = "|"))) %>% select(c(contains("query"),"strand","gene")) %>% mutate(gene = gene %>% str_remove_all(".[0-9].v3.1"))%>%
    left_join(filter(subgenome_key, Genome == tripsacinae_genome_IDs[t]), by="queryChr")
  QC.genes[,paste0(tripsacinae_genome_IDs[t],"_chr")]<-NA
  QC.genes[,paste0(tripsacinae_genome_IDs[t],"_start")]<-NA
  QC.genes[,paste0(tripsacinae_genome_IDs[t],"_end")]<-NA
  QC.genes[,paste0(tripsacinae_genome_IDs[t],"_strand")]<-NA
  for(i in 1:nrow(QC.genes)){
    parsing<-filter(subgenome_key, Genome == tripsacinae_genome_IDs[t] & RefChr == QC.genes$CHROM[i] & Subgenome == QC.genes$Subgenome[i])
    if(nrow(filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i]& RefChr == QC.genes$CHROM[i])) >0){
      QC.genes[i,paste0(tripsacinae_genome_IDs[t],"_chr")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,1]
      QC.genes[i,paste0(tripsacinae_genome_IDs[t],"_start")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,2]
      QC.genes[i,paste0(tripsacinae_genome_IDs[t],"_end")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,3]
      QC.genes[i,paste0(tripsacinae_genome_IDs[t],"_strand")]<-filter(df, gene == QC.genes$Gene_ID[i] & queryChr %in% parsing$queryChr & Subgenome == QC.genes$Subgenome[i] & RefChr == QC.genes$CHROM[i])[,4]
    }
  }
  print(paste("Finished with",t,"/35 genomes :)"))
}

#Pick out genes that are used as anchors across most genomes (30/35?)
QC.genes<-mutate(QC.genes, CntNA = is.na(pick(contains("_start"))) %>% rowSums()) #adds a column with the count of genomes that are NA for that gene/subgenome

QC.genes.clean<-filter(QC.genes, CntNA <= 5) #Goes from 24320 to 10665
#QC.genes.clean<-filter(QC.genes, CntNA <= 15) #Goes from 24320 to 1452

filter(QC.genes, CHROM %in% c("Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")) %>% group_by(CHROM) %>% summarize(minNA = min(CntNA), maxNA = max(CntNA))

#Pick out a few random genes from the genes with mixed alignment calls
QC.Gene_IDs<-mutate(genes_with_mixed_alignment, Gene_ID = str_remove(string=genes_with_mixed_alignment$Gene_ID, pattern = ".[0-9].v3.[0-9]")) %>%filter(Gene_ID %in% QC.genes.clean$Gene_ID) %>% select(Gene_ID) %>% unique() %>% pull()%>% sample(.,4) 
#There are none that pass the filtering above that are also in the mixed alignment

#Pick out a few random genes from each ref chr
for(i in c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10")){
  if(nrow(filter(QC.genes.clean, CHROM == i)) > 0 ){
    QC.Gene_IDs<-c(QC.Gene_IDs, 
                   filter(QC.genes.clean, CHROM == i) %>% select(Gene_ID) %>% pull() %>% sample(.,4))
  }
}
QC.Gene_IDs<-unique(QC.Gene_IDs)
#check that length == 50 for unique gene IDs
length(QC.Gene_IDs)
#Rerun until length is at least 50 or run below
for(i in c("Chr10")){
  if(nrow(filter(QC.genes.clean, CHROM == i)) > 0 ){
    QC.Gene_IDs<-c(QC.Gene_IDs, 
                   filter(QC.genes.clean, CHROM == i) %>% select(Gene_ID) %>% pull() %>% sample(.,10))
  }
}
QC.Gene_IDs<-unique(QC.Gene_IDs)
length(QC.Gene_IDs)
#If more than 50, head the first 50 lines
QC.Gene_IDs<-head(QC.Gene_IDs, n=50)
QC.Gene_IDs<-str_remove_all(QC.Gene_IDs, ".[0-9].v3.1")

#Pull out the subset Gene IDs to be a table to send to Maggie
#Table 1 coordinates
write_tsv(filter(QC.genes.clean, Gene_ID %in% QC.Gene_IDs), "/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/QC/QC.Gene_IDs.Coordinates.tsv")
mutate(full.fractionation.status, Gene_ID = str_remove(string=full.fractionation.status$Gene_ID, pattern = ".[0-9].v3.[0-9]")) %>%filter(Gene_ID %in% QC.Gene_IDs) %>%
  select(c(contains("ID"),ends_with(".status"))) %>% write_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/QC/QC.Gene_IDs.FractionationStatus.tsv")

#Pull out deletion lengths for deletions that intersect the QC genes
#In bash
#ml bedtools2
# grep "Chr10" QC.Gene_IDs.Coordinates.tsv | awk -v OFS='\t' '{print $1,$2,$3,$7,$5}'| sort -k2n,2n | sed 's/Chr//g' | bedtools intersect -wb -a - -b ../tripsacinae-sb_split_mafs/Sb313_TdFL_Chr10.bin1.dels.bed | head -n 1 |awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$7,$8,$13}'
# echo Sbgene_Chr Sbgene_Start Sbgene_End Sbgene_Strand Gene_ID Del_Start Del_Stop Del_Length Genome Subgenome | tr " " "\t" > QC.deletions.tsv
#for g in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML333 ZmCML322 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
#for c in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 ; do 
#grep $c QC.Gene_IDs.Coordinates.tsv | awk -v OFS='\t' '{print $1, $2, $3, $7, $5}' | sort -k2n,2n | sed 's/Chr//g' | bedtools intersect -wb -a - -b ../tripsacinae-sb_split_mafs/Sb313_${g}_${c}.bin1.dels.bed | awk -v OFS='\t' -v genome=$(echo $g) '{print $1,$2,$3,$4,$5,$7,$8,$13,genome, "M1"}' >> QC.deletions.tsv ; 
#grep $c QC.Gene_IDs.Coordinates.tsv | awk -v OFS='\t' '{print $1, $2, $3, $7, $5}' | sort -k2n,2n | sed 's/Chr//g' | bedtools intersect -wb -a - -b ../tripsacinae-sb_split_mafs/Sb313_${g}_${c}.bin2.dels.bed | awk -v OFS='\t' -v genome=$(echo $g) '{print $1,$2,$3,$4,$5,$7,$8,$13,genome, "M2"}' >> QC.deletions.tsv ; 
#done ; echo Done with ${g} ; done

#What is the distribution of deletion lengths across genomes?
for(i in tripsacinae_genome_IDs){
  for(c in c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")){
    for(m in c("bin1","bin2")){
      assign(paste0("temp_",i,"_",c,"_",m), 
             read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/Sb313_",i,"_",c,".",m,".dels.bed"), 
                      col_names = c("Chr","Start","End","REF","ALT","Strand","SingularALT", "DEL_length", "genotype")))
    }
  }
}

#Use rbind() to append dataframes
#Test before running on everything
for(i in tripsacinae_genome_IDs){
  for(b in c("bin1","bin2")){
    assign(paste0(i,"_",b), get(paste0("temp_",i,"_Chr01_",b)))
    for(c in c("Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10")){
      assign(paste0(i,"_",b), rbind(get(paste0(i,"_",b)),get(paste0("temp_",i,"_",c,"_",b))))
    }
  }
}

remove(list = ls(pattern = "temp_")) #to free up some space
#add genome and subgenome columns to each dataframe
for(i in tripsacinae_genome_IDs){
  for(b in c("bin1","bin2")){
    assign(paste0(i,"_",b), mutate(get(paste0(i,"_",b)), Genome = i, Subgenome = case_when(b == "bin1" ~ "M1", b == "bin2" ~ "M2")))
  }
}
#combine all dataframes (and thus all deletions called on every genome) into a single dataframe
for(i in tripsacinae_genome_IDs){
  for(b in c("bin1","bin2")){
    if(i == tripsacinae_genome_IDs[1] & b == "bin1"){
      all_deletions<-get(paste0(i,"_",b))
    }else{
      all_deletions<-rbind(all_deletions, get(paste0(i,"_",b)))
    }
  }
}
#make space by removing the extra dataframes
remove(list = ls(pattern = "_bin1")) #to free up some space
remove(list = ls(pattern = "_bin2"))

##Plot distribution of deletion lengths
#deletions by subgenome
ggplot(all_deletions, aes(x=Subgenome, y=DEL_length))+
  geom_boxplot(aes(fill = Subgenome))+
  scale_fill_manual(values = subgenome_colors)+
  theme_bw()+
  ylab("Deletion Length (bp)")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/alldelLengths.bySubgenome.png", device="png",dpi=300, height=6, width=6)

#deletions by genome
ggplot(all_deletions, aes(x=Genome, y=DEL_length))+
  geom_boxplot(aes(fill = Genome))+
  scale_fill_manual(values = genome_colors)+
  theme_bw()+
  ylab("Deletion Length (bp)")+
  theme(axis.text.x = element_text(angle = 90))
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/alldelLengths.byGenome.png", device="png",dpi=300, height=6, width=)


#deletions by subgenome and genome
ggplot(all_deletions, aes(x=Subgenome, y=DEL_length, fill=Genome))+
  geom_boxplot()+
  scale_fill_manual(values=genome_colors)+
  theme_bw()+
  ylab("Deletion Length (bp)")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/alldelLengths.bySubgenomeAndGenome.png", device="png",dpi=300, height=6, width=6)


na.omit(all_deletions) %>% group_by(Genome,Subgenome) %>% reframe(min.del = min(DEL_length),
                                                         max.del = max(DEL_length),
                                                         mean.del = mean(DEL_length),
                                                         median.del = median(DEL_length))