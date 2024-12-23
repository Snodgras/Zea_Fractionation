#TROUBLESHOOTING
library(GenomicRanges)
library(tidyverse)

#Step 0: Filter by ref chromosome, query genome name, and query x ref blocks; make into gr object
step0<-function(REFchr,queryGENOME,refGENOME){ #REFchr = "Chr10", queryGENOME = "ZmB73",refGENOME = "Sb313"
  uniqblkID<-filter(phased_blk, refChr == REFchr & grepl(queryGENOME, blkID)) %>% select(blkID) %>% unique() #create a vector with their unique BlockID strings
  df<-filter(phased_blk, refChr == REFchr & grepl(queryGENOME, blkID)) %>% arrange(.$blkID) #make it into a dataframe
  df<-df[match(uniqblkID$blkID,df$blkID),]
  df<-filter(df, genome1 == refGENOME | genome2 == refGENOME)
  return(df)
}

#t1<-step0("Chr07","ZmB73","Sb313") #I will supply this
t1<-read_tsv("Chr07.parsingdata.tsv")

#Step 1: Split gr object by query chromosomes
step1<-function(df){
  gr<-df %>%
    GRanges(seqnames = Rle(.$chr1),
            ranges = IRanges(start = .$startBp1, end = .$endBp1, names = .$blkID),
            strand = Rle(strand(.$orient)),
            queryChr = .$chr2,
            queryStart = .$startBp2,
            queryEnd = .$endBp2)
  s.gr<-split(gr, gr$chr2) #splits the genomic range object by the query chr
  return(s.gr)
}

s.t1<-step1(t1)

#Step 2: Find the query chromsome that contains the most ref bps (max width)
step2<-function(s.df){
  t<-c()
  for(i in 1:length(s.df)){
    t<-c(t,sum(width(ranges(s.df)[i]) %>% unlist()))
  }
  #need to find the max length
  return(grep(max(t),t))
}

y<-step2(s.t1)

#Step 3: Overlap the ranges of the smaller width chromosomes with the widest chr from 2
#save overlap as a gr object
#save the widest chr as gr object
#if more than 3 query chromosomes, check to make sure there's no overlap between the smaller chromosome chunks

step3.a<-function(s.df, y){ #checks for overlaps between large chunk (y) and smaller chunks (i)
  for(i in 1:length(s.df)){
    if(i != y){
      assign(paste("int",i,"by",y, sep = "."),subsetByOverlaps(s.df[i],s.df[y], ignore.strand = TRUE) %>% unlist(), envir = parent.frame())
    }
    else{
      assign("dup.1", unlist(s.df[y]), envir = parent.frame())
    }
  }
}
step3.b<-function(s.df, y){
  if(length(s.df) > 2){ #first check that this is even needed by checking if there are more than 2 query chromosomes for this region
    l<-c(1:length(s.df))
    l<-l[!grepl(y,l)]
    if(length(l) == 2){
      if(subsetByOverlaps(s.df[l[1]],s.df[l[2]],ignore.strand = TRUE) %>% length() != 0){
        assign(paste("int",l[1],"by",l[2], sep = "."), subsetByOverlaps(s.df[l[1]],s.df[l[1]],ignore.strand = TRUE) %>% unlist(),envir = parent.frame())
        print("STEP3b: There are overlaps between the smaller chunks; may need clean up")
      }
      else{print("STEP3b: There are no overlaps between the smaller chunks; this is good :)")}
    }
    else{print("STEP3b: There are more than 2 small chunks :|")
      for(a in 1:length(l)){
        for(b in 1:length(l)){
          if(subsetByOverlaps(s.df[l[a]],s.df[l[b]],ignore.strand = TRUE) %>% length() != 0 & a != b ){
            assign(paste("int",l[a],"by",l[b], sep = "."), subsetByOverlaps(s.df[l[a]],s.df[l[b]],ignore.strand = TRUE) %>% unlist(),envir = parent.frame())
            print("STEP3b: There are overlaps between the smaller chunks; may need clean up")
          } 
        }
      }
    }
  }else{print("STEP3b: There is only <= 1 small chunk, no need to check for overlaps :)")}
}
step3.a(s.t1,y) #checks for large chunk overlap with smaller chunks
step3.b(s.t1,y) #checks for overlap between the smaller chunks

#Step 4: Save the gr objects from 3 as temp dataframes with specific column names to easily make a bed file
# widest chromosomes as file 2
# overlapped chrs as file 2
step4.a<-function(int.df){ #makes a single dataframe
  dup.regions<-data.frame(seqnames=seqnames(int.df),
                          starts=start(int.df)-1, #-1 because bed files are 0 indexed and GR is 1 indexed
                          ends=end(int.df),
                          names=int.df$blkID,
                          scores=rep(".", length(int.df)),
                          strands=strand(int.df),
                          query_chr=int.df$queryChr,
                          query_start=int.df$queryStart-1,
                          query_end=int.df$queryEnd)
  return(dup.regions)
}
step4.b<-function(df1, df2){ #combines dataframes if needed
  dup.regions<-add_row(df1, seqnames = df2$seqnames,
                       starts = df2$starts,
                       ends = df2$ends,
                       names = df2$names,
                       scores = df2$scores,
                       strands = df2$strands,
                       query_chr = df2$query_chr,
                       query_start = df2$query_start,
                       query_end = df2$query_end)
  return(dup.regions)
}

#I made the use of this function into a loop so I don't have to manually specify file names
for(k in ls(pattern=".by.")){
  assign(paste0(k,".df",sep=""), step4.a(get(k))) #convert to df each intersect file
}
dup1.df<-step4.a(dup.1)

####PROBLEM IS ANYTHING THAT IS NOT THE MAIN CHROMSOME PIECE IS PUT INTO DUP2, EVEN IF THERE'S OVERLAP IN THE SB313 COORDS

if(length(ls(pattern = ".by.[0-9].df")) > 1){ #if there's multiple intersect dfs, combine them
  for(m in ls(pattern = ".by.[0-9].df")){
    if(exists("dup.2.df")){
      dup.2.df<-step4.b(dup.2.df,get(m)) #if the dataframe already exists, add to it
    }
    else{dup.2.df<-get(m)} #if it doesn't exist, make it with the first int.df
  }
}else{dup.2.df <-get(ls(pattern = ".by.[0-9].df"))} #if there's not multiple intersect dfs

#Step 5: Quality control the output to make sure there's no errors
#QC checks:
#1. Total n rows of the dup.regions dataframes == nrow of the filtered df
#2. All blockIDs found in the filtered df are also found in the duplicated regions:
#3. All blockIDs found in the duplicated region files are in there exactly once

step5.a<-function(orig.df,dup1,dup2){
  QC1<-nrow(orig.df) == sum(c(nrow(dup1),nrow(dup2)))
  return(QC1)
}
step5.b<-function(orig.df,dup1,dup2){
  QC2<-filter(orig.df) %>% select("blkID") %>%
    filter(!blkID %in% c(dup1$names,dup2$names)) %>% nrow() == 0
  return(QC2)
}
step5.c<-function(orig.df,dup1,dup2){
  QC3<-c(duplicated(dup1$names) %>% sum(), #should sum to 0 if all are FALSE (means no duplicates)
         duplicated(dup2$names) %>% sum(),
         duplicated(c(dup1$names,dup2$names))) %>% sum()
  return(QC3)
}
step5.a(t1,dup.1.df,dup.2.df) #TRUE TO PASS
step5.b(t1,dup.1.df,dup.2.df) #TRUE TO PASS
step5.c(t1,dup.1.df,dup.2.df) #0 TO PASS

#Step 6: Add the temp dataframes to a final dataframe
step6<-function(final.df,dup.df){
  if(exists(final.df)){
    df<-add_row(get(final.df), seqnames = dup.df$seqnames,
                starts = dup.df$starts,
                ends = dup.df$ends,
                names = dup.df$names,
                scores = dup.df$scores,
                strands = dup.df$strands,
                query_chr = dup.df$query_chr,
                query_start = dup.df$query_start,
                query_end = dup.df$query_end)
  }else{
    df<-dup.df
  }
  return(df)
}
Sb313_vs_ZmB73.dupregions.1<-step6("Sb313_vs_ZmB73.dupregions.1",dup.1.df)
Sb313_vs_ZmB73.dupregions.2<-step6("Sb313_vs_ZmB73.dupregions.2",dup.2.df)
