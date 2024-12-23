library(tidyverse)

# Function to process convergence for a single gene
runConvergencePerGene <- function(gene_id, output_dir) {
  ref.gene <- read_lines("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/ref.gene.txt")
  tripsacinae_genome_IDs <- read_lines("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/tripsacinae_genome_IDs.txt")
  gene_atleastOneNA.NoZnZd <- read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/gene_atleastOneNA.NoZnZd.tsv")
  convergence.NoZnZd <- read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/convergence.NoZnZd.tsv")

  convergence_shared <- tibble(Gene_ID = NA, Subgenome = NA,
                               Target_Genome = NA, Query_Genome = NA,
                               Convergence_Category = NA)

  for (m in c("M1", "M2")) { # for each subgenome
    if (nrow(filter(gene_atleastOneNA.NoZnZd, M == m & Gene_ID == gene_id)) == 0) {
      df <- filter(convergence.NoZnZd, Gene_ID == gene_id & M == m)
      
      for (g in tripsacinae_genome_IDs[c(1:4, 8:24, 27:32, 34, 36:39)]) {
        pattern1 <- filter(df, Genome == g) %>% select(Loss_Pattern) %>% pull()

        for (q in tripsacinae_genome_IDs[c(1:4, 8:24, 27:32, 34, 36:39)]) {
          pattern2 <- filter(df, Genome == q) %>% select(Loss_Pattern) %>% pull()

          if (is_empty(c(pattern1, pattern2))) {
            convergence_shared <- add_row(convergence_shared, Gene_ID = gene_id, Subgenome = m,
                                          Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyShared_Retention")
          } else if ((is_empty(pattern1) | is_empty(pattern2)) & !is_empty(c(pattern1, pattern2))) {
            convergence_shared <- add_row(convergence_shared, Gene_ID = gene_id, Subgenome = m,
                                          Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyDifferent_Retention")
          } else if (pattern1 == pattern2) {
            convergence_shared <- add_row(convergence_shared, Gene_ID = gene_id, Subgenome = m,
                                          Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyShared")
          } else {
            if (length(str_split(pattern1, ":", simplify = T)) >= length(str_split(pattern2, ":", simplify = T))) {
              longest.pattern <- pattern1
              shortest.pattern <- pattern2
            } else {
              longest.pattern <- pattern2
              shortest.pattern <- pattern1
            }

            shortest.pattern <- str_split(shortest.pattern, ":", simplify = T)
            counter <- 0
            for (e in 1:length(shortest.pattern)) {
              if (str_detect(longest.pattern, shortest.pattern[1, e])) {
                counter <- counter + 1
              }
            }

            if (counter > 0) {
              convergence_shared <- add_row(convergence_shared, Gene_ID = gene_id, Subgenome = m,
                                            Target_Genome = g, Query_Genome = q, Convergence_Category = "SomeShared")
            } else {
              convergence_shared <- add_row(convergence_shared, Gene_ID = gene_id, Subgenome = m,
                                            Target_Genome = g, Query_Genome = q, Convergence_Category = "CompletelyDifferent")
            }
          }
        }
      }
    }
  }

  output_file <- file.path(output_dir, paste0(gene_id, "_convergence_shared.tsv"))
  write_tsv(convergence_shared, output_file)
}

