library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggtranscript) #will need to add to container
library(GenomicRanges)
library(scales)
#setwd(
#  "C:/Users/arun/OneDrive - Iowa State University/OrganizedDocuments/HuffordLab/Sam/bedfile"
#)
Sb313_phasedBlks <- read.csv("/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/Sb313_phasedBlks.csv", header=TRUE)

# for Sorghum
sb_filtered <- Sb313_phasedBlks %>%
  filter(genome1 == "Sb313") %>%
  mutate(segment_name = paste0(chr2, ":", startBp2, "-", endBp2)) %>%
  select(genome2, 
         chr1, 
         startBp1, 
         endBp1, 
         segment_name,
         orient) %>%
  filter(str_detect(chr1, "^Chr")) %>%
  mutate(group = gsub('(.*):.*', '\\1', segment_name))

# rename header
colnames(sb_filtered) <- c("query_genome",
                           "seqnames",
                           "start",
                           "end",
                           "segment_name",
                           "strand", 
                           "group")
# # for B73
# b73_filtered <- Sb313_phasedBlks %>%
#   filter(genome1 == "ZmB73") %>%
#   mutate(segment_name = paste0(chr2, ":", startBp2, "-", endBp2)) %>%
#   select(genome2, 
#          chr1, 
#          startBp1, 
#          endBp1, 
#          segment_name,
#          orient) %>%
#   filter(str_detect(chr1, "^chr"))  %>%
#   mutate(group = gsub('(.*):.*', '\\1', segment_name))
# # rename header
# colnames(b73_filtered) <- c("query_genome",
#                            "seqnames",
#                            "start",
#                            "end",
#                            "segment_name",
#                            "strand", 
#                            "group")
# separate each spp alignment and remove the query_genome col
# for sorghum
sb_filteredList <-
  lapply(split(sb_filtered, sb_filtered$query_genome), function(x) {
    x$query_genome <- NULL
    x
  })
# for B73
# b73_filteredList <-
#   lapply(split(b73_filtered, b73_filtered$query_genome), function(x) {
#     x$query_genome <- NULL
#     x
#   })
# # B73 length file for GRanges   
# b73_length <- read.table("b73_length.tsv", header= FALSE, sep = "\t")
# b73_chrLen <- Seqinfo(seqnames = b73_length$V1,
#                      isCircular = rep(FALSE, length(b73_length$V1)),
#                      seqlengths=b73_length$V2,
#                      genome="zmB73v5")
# sorghum length file for GRanges   
sb_length <- read.table("/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/sb_length.bed", header= FALSE, sep = "\t")
sb_chrLen <- Seqinfo(seqnames = sb_length$V1,
                     isCircular = rep(FALSE, length(sb_length$V1)),
                     seqlengths=sb_length$V2,
                     genome="sb313")

# gr.names.refB73 <- vector("character")
# for (i in 1:length(b73_filteredList)) {
#   varName <- names(b73_filteredList)[i]
#   myName <- paste0("zmB73vs", varName)
#   gr.names.refB73[length(gr.names.refB73) + 1] = myName
#   myGR <- makeGRangesFromDataFrame(
#     b73_filteredList[[i]],
#     keep.extra.columns = TRUE,
#     ignore.strand = FALSE,
#     seqinfo = b73_chrLen,
#     seqnames.field = "seqnames",
#     start.field = "start",
#     end.field = c("end", "stop"),
#     strand.field = "strand",
#     starts.in.df.are.0based = TRUE
#   )
#   gList = split(myGR, seqnames(myGR))
#   assign(myName, gList)
# }

gr.names.refSB <- vector("character")
for (i in 1:length(sb_filteredList)) {
  varName <- names(sb_filteredList)[i]
  myName <- paste0("sb313vs", varName)
  gr.names.refSB[length(gr.names.refSB) + 1] = myName
  myGR <- makeGRangesFromDataFrame(
    sb_filteredList[[i]],
    keep.extra.columns = TRUE,
    ignore.strand = FALSE,
    seqinfo = sb_chrLen,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = c("end", "stop"),
    strand.field = "strand",
    starts.in.df.are.0based = TRUE
  )
  gList = split(myGR, seqnames(myGR))
  assign(myName, gList)
}

mergeGRobj <- function(grObj = gChr10.split, i = 1) {
  grns <-
    GRanges(seqnames = seqnames(grObj[[i]]), 
            ranges = ranges(grObj[[i]]))
  name <- names(grObj)[i]
  outGR <- reduce(grns, min.gapwidth = 50000000)
  outGR$group <- name
  return(outGR)
}

# delete sb313vsAv and sb313vsSb313
gr.names.refSB <- gr.names.refSB[3:length(gr.names.refSB)]

for (grName in gr.names.refSB) {
  gList = get(grName)
  qName <- str_sub(grName, 8,)
  results <- paste0("binned.", grName)
  for (name in names(gList)) {
    item <- paste0("gr.", name)
    item.merged <- paste0("gr.", name, ".merged")
    item.bin1 <- paste0("gr.", name, ".bin1")
    item.bin2 <- paste0("gr.", name, ".bin2")
    split.genome <- paste0(name, ".", grName)
    gr.split.merged = list()
    gr.split <- split(gList[[name]], gList[[name]]$group)
    for (i in 1:length(gr.split)) {
      output <- mergeGRobj(grObj = gr.split, i = i)
      gr.split.merged[[length(gr.split.merged) + 1]] = output
    }
    names(gr.split.merged) <-
      sapply(gr.split.merged, function(gr.split)
        gr.split$group)
    largest_segment <- NULL
    largest_length <- 0
    largest_name <- NULL
    other_segments <- list()
    for (i in 1:length(gr.split.merged)) {
      current_length <- width(gr.split.merged[[i]])
      currnet_name <- gr.split.merged[[i]]$group
      if (current_length > largest_length) {
        largest_length <- current_length
        largest_name <- currnet_name
        largest_segment <- gr.split.merged[[i]]
      }
    }
    other_segments <-
      gr.split.merged[!sapply(gr.split.merged, identical, largest_segment)]
    bin1 <- largest_segment$group
    bin2 <- NULL
    for (i in 1:length(other_segments)) {
      overlaps <- findOverlaps(largest_segment, other_segments[[i]])
      if (subjectHits(overlaps) > 0) {
        bin2[length(bin2) + 1] <- other_segments[[i]]$group
      } else {
        bin1[length(bin1) + 1] <- other_segments[[i]]$group
      }
    }
    gr.split.bin1 <- unlist(gr.split[bin1])
    gr.split.bin2 <- unlist(gr.split[bin2])
    gr.split.bin1.df <-
      as.data.frame(gr.split.bin1, row.names = NULL)
    gr.split.bin2.df <-
      as.data.frame(gr.split.bin2, row.names = NULL)
    gr.SB.bin1 <- gr.split.bin1.df %>%
      separate(segment_name,
               into = c('qryChr', 'qryTemp'),
               sep = ':') %>%
      separate(qryTemp,
               into = c('qryStart', 'qryEnd'),
               sep = '-') %>%
      mutate(ref.segment = paste0(seqnames, ":", start, "-", end)) %>%
      mutate(ref.group = seqnames) %>%
      select(qryChr, qryStart, qryEnd, strand, ref.segment, ref.group)
    
    gr.SB.bin2 <- gr.split.bin2.df %>%
      separate(segment_name,
               into = c('qryChr', 'qryTemp'),
               sep = ':') %>%
      separate(qryTemp,
               into = c('qryStart', 'qryEnd'),
               sep = '-') %>%
      mutate(ref.segment = paste0(seqnames, ":", start, "-", end)) %>%
      mutate(ref.group = seqnames) %>%
      select(qryChr, qryStart, qryEnd, strand, ref.segment, ref.group)
    
    
    gr.SB.bin2 <- as_tibble(gr.SB.bin2)
    gr.SB.bin2$qryStart <- as.numeric(gr.SB.bin2$qryStart) 
    gr.SB.bin2$qryEnd <- as.numeric(gr.SB.bin2$qryEnd) 
    
    gr.SB.bin1 <- as_tibble(gr.SB.bin1)
    gr.SB.bin1$qryStart <- as.numeric(gr.SB.bin1$qryStart) 
    gr.SB.bin1$qryEnd <- as.numeric(gr.SB.bin1$qryEnd) 
    
    gr.SB.bin1$bin <- "bin1"
    gr.SB.bin2$bin <- "bin2"
    Sb_Qry <- rbind(gr.SB.bin1, gr.SB.bin2)
    assign(split.genome, Sb_Qry)
  }
  regex_pattern <- paste0("^Chr.*", grName)
  object_names <- ls(pattern = regex_pattern)
  object_list <- list()
  for (name in object_names) {
    object_list[[name]] <- get(name)
  }
  assign(results, object_list)
}

