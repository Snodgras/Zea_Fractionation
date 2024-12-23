library(tidyverse)

#Find genes for gene expression
#Genes that have different exons fractionated between NAM lines

gene_expression_subset<-filter(convergence.NoZnZd, str_detect(Genome, "Zm"))

genes_multLossPatterns<-gene_expression_subset %>% 
  group_by(Gene_ID, M) %>% 
  select(Gene_ID, M, Loss_Pattern) %>% 
  unique() %>% 
  count() %>% 
  filter(n > 1)

#just to see how many loss patterns we'd be comparing for a given gene
#looks to be 2-16 with the majority being 2-3
summary(genes_multLossPatterns)

gene_expression_subset<-filter(gene_expression_subset, Gene_ID %in% genes_multLossPatterns$Gene_ID )

#Creates table/df where each row is homoeolog and each genome a column with the loss pattern as the value
gene_expression_subset_pattern_table<-gene_expression_subset %>% select(Gene_ID, ExonCt, Genome, M, Loss_Pattern) %>% 
  pivot_wider(names_from = "Genome", values_from = "Loss_Pattern")

#Create M1 and M2 lists of Gene_IDs
filter(gene_expression_subset, M == "M1") %>% select(Gene_ID) %>% pull() %>% unique() %>% str_remove_all("\\.v3\\.[0-9]*") %>% write_lines("/work/LAS/mhufford-lab/snodgras/Fractionation/gene_expression/M1_Sb313_geneIDs.txt")
filter(gene_expression_subset, M == "M2") %>% select(Gene_ID) %>% pull() %>% unique() %>% str_remove_all("\\.v3\\.[0-9]*") %>% write_lines("/work/LAS/mhufford-lab/snodgras/Fractionation/gene_expression/M2_Sb313_geneIDs.txt")

# Sorghum - NAM gene names from https://ars-usda.app.box.com/v/maizegdb-public/folder/186350887665/Sorghum_NAM_synteny_associations_subgenomes_longest_transcripts_bed_files.tar.gz
# NAM TPM gene expression from https://ars-usda.app.box.com/v/maizegdb-public/folder/175926197781?page=1

#In unix
## Get the gene IDs for each unique genome
# for genome in Sorghum_NAM_synteny_associations_subgenomes_longest_transcripts_bed_files/*.bed ; do 
#   gnam=$(echo ${genome} | cut -f 2 -d "-")
#   echo ${gnam}
#   grep "M1" ${genome} | grep -f M1_Sb313_geneIDs.txt | cut -f 4-5 > ${gnam}_M1_gene_IDs.tsv
#   grep "M2" ${genome} | grep -f M2_Sb313_geneIDs.txt | cut -f 4-5 > ${gnam}_M2_gene_IDs.tsv
# done

## Get the expression levels for each unique genome for subsets of gene IDs
# cd NAM_pangene_expression_counts_per_tissue-TPM/
# mv Oh7b_full-tpm.tsv Oh7B_full-tpm.tsv
# mv IL14H_full-tpm.tsv Il14H_full-tpm.tsv
# mv MS71_full-tpm.tsv Ms71_full-tpm.tsv
# for gene in *M*gene_IDs.tsv ; do 
# gnam=$(echo ${gene} | cut -f 1 -d "_")
# subgen=$(echo ${gene} | cut -f 2 -d "_")
# echo ${gnam} ${subgen}
# head -n 1 NAM_pangene_expression_counts_per_tissue-TPM/${gnam}_full-tpm.tsv > ${gnam}_${subgen}_subset-tpm.tsv
# cut -f 1 ${gene} | grep -f - NAM_pangene_expression_counts_per_tissue-TPM/${gnam}_full-tpm.tsv >> ${gnam}_${subgen}_subset-tpm.tsv
#done

# Create Key of gene IDs
NAM_ids<-c("B73","B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301","Il14H","Ki11","Ki3","Ky21","M162W","M37W","Mo18W","Ms71","NC350","NC358","Oh43","Oh7B","P39","Tx303","Tzi8")

for(i in NAM_ids){
  assign(paste0(i,"_Sb_M1_geneIDs"), read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/gene_expression/",i,"_M1_gene_IDs.tsv"), col_names = c(paste0(i,"_GeneID"),"Sb313_GeneID")))
  assign(paste0(i,"_Sb_M2_geneIDs"), read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/gene_expression/",i,"_M2_gene_IDs.tsv"), col_names = c(paste0(i,"_GeneID"),"Sb313_GeneID")))
}

#bring in expression data
for(i in NAM_ids){
  assign(paste0(i,"_M1_TPM"), read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/gene_expression/",i,"_M1_subset-tpm.tsv")))
  assign(paste0(i,"_M2_TPM"), read_tsv(paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/gene_expression/",i,"_M2_subset-tpm.tsv")))
}

#left join to get sorghum ids into the expression dataframes
left_join(x=B73_M1_TPM, y=B73_Sb_M1_geneIDs, by=c("Geneid"="B73_GeneID"))
for(i in NAM_ids){
  assign(paste0(i,"_M1_TPM"), left_join(x=get(paste0(i,"_M1_TPM")), y=get(paste0(i,"_Sb_M1_geneIDs")), by=c("Geneid" = paste0(i,"_GeneID"))))
  assign(paste0(i,"_M2_TPM"), left_join(x=get(paste0(i,"_M2_TPM")), y=get(paste0(i,"_Sb_M2_geneIDs")), by=c("Geneid" = paste0(i,"_GeneID"))))
}

#Compare across genomes using tissues+individuals as replicates
gene_expression_subset_pattern_table[1,] %>% 
  pivot_longer(cols = starts_with("Zm"), names_to = "Genome", values_to = "Pattern") %>%
  print(n=26)

gene_expression_subset_pattern_table %>% 
  pivot_longer(cols = starts_with("Zm"), names_to = "Genome", values_to = "Pattern") %>%
  group_by(Gene_ID, M, Pattern) %>% 
  count() %>% 
  filter(n!=26) %>%
  ggplot(aes(x=n))+geom_histogram(binwidth = 1)
