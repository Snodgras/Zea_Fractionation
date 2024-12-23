#!/usr/bin/env Rscript

#Example script that needs to be parallelized
library(tidyverse)

ref.gene<-read_lines("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/ref.gene.txt")
tripsacinae_genome_IDs<-read_lines("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/tripsacinae_genome_IDs.txt")
gene_atleastOneNA.NoZnZd<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/gene_atleastOneNA.NoZnZd.tsv")
convergence.NoZnZd<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/convergence.NoZnZd.tsv")

convergence_shared<-tibble(Gene_ID =NA, Subgenome = NA, 
                           Target_Genome = NA, Query_Genome = NA, 
                           Convergence_Category = NA)

for(i in ref.gene){#for each gene
  for(m in c("M1","M2")){ #for each subgenome
    #so long as the gene x subgenome combo isn't in the at least one NA dataframe
    #i.e. all genomes included have that gene as aligned
    if(nrow(filter(gene_atleastOneNA.NoZnZd, M == m & Gene_ID == i)) == 0){
      #filter the larger convergence data frame to just the gene x subgenome we want to work with
      df<-filter(convergence.NoZnZd, Gene_ID == i & M == m)
      #for each genome as target
      for(g in tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]){
        #pull out the loss pattern for target genome
        pattern1 <- filter(df, Genome == g) %>% select(Loss_Pattern) %>% pull()
        #now for each query genome (including the self comparison) 
        for(q in tripsacinae_genome_IDs[c(1:4,8:24,27:32,34,36:39)]){
          #pull out the loss pattern for the query genome
          pattern2 <- filter(df, Genome == q) %>% select(Loss_Pattern) %>% pull()
          #if both patterns are empty (meaning no exon was called fractionated for either genome for gene)
          if(is_empty(c(pattern1,pattern2))){
            #call that completely shared (by) Retention
            convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                        Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyShared_Retention")
          }else{
            #if one of the patterns is empty but not both
            if((is_empty(pattern1) | is_empty(pattern2)) & !is_empty(c(pattern1,pattern2))){
              #call that completely different (by) Retention
              convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                          Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyDifferent_Retention")
              
            }else{
              #if neither pattern is empty
              if(!is_empty(c(pattern1,pattern2))){
                #then check for matching of pattern
                if(pattern1 == pattern2){ #if the patterns match exactly, call that Completely Shared
                  convergence_shared<-add_row(convergence_shared,Gene_ID =i, Subgenome = m, 
                                              Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyShared")
                }else{
                  #figure out which is the longest pattern, but only if both patterns aren't empty (extra precaution given previous if/else statements)
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

write_tsv(convergence_shared, "/work/LAS/mhufford-lab/snodgras/Fractionation/intermediate-data-files/convergence_shared.tsv")
