library(GENESPACE)
#spec <- c("zTIL01", "zTIL11", "zTIL18", "zTIL25", "zB73v5", "zdgigi", "zdmomo", "zmhuet", "avirgi", "ppanic", "sbicol", "rtuber", "tdacts", "irugos" )
#ploidy <- c(2,2,2,2,2,2,2,2,1,1,1,1,2,1)
#ref <- "sbicol"

#names is a tsv file with the first field the species name (abbrev.) and second field is ploidy
#WILL WE NEED TO SPECIFY THE SORGHUM AS DIPLOID AND EVERYTHING ELSE AS TETRAPLOID (YES)
#for prelim run
#names <- read.delim("/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/names.tsv")
#reqNames <- c("Sb313","TdFL","ZdGigi","ZvTIL01", "ZxTIL18","ZmB73","ZmCML333","ZmNC358","ZmOh43","Av")
#for full run
names <- read.delim("/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2/names.tsv")
reqNames <- c("Sb313", "TdFL","TdKS","ZdGigi","ZdMomo","ZhRIMHU001",
              "ZnPI615697",
              "ZvTIL01","ZvTIL11",
              "ZxTIL18","ZxTIL25","ZmB73","ZmB97","ZmCML103","ZmCML228","ZmCML247","ZmCML277","ZmCML322",
              "ZmCML333","ZmCML52","ZmCML69","ZmHP301","ZmIL14H",
              "ZmKi11","ZmKi3","ZmKy21","ZmM162W",
              "ZmM37W","ZmMo18W","ZmMS71",
              "ZmNC350","ZmNC358","ZmOh43","ZmOh7b",
              "ZmP39","ZmTx303",
              "ZmTzi8","Av")

#Specifies the genomes we want to plot so that it speeds up plotting
#Can specify anything we want
#To run full pairs, then we'd use the full table in the namesdf
namesdf <-  names[names$spec %in% reqNames ,]

spec <- namesdf$spec
ploidy <- namesdf$ploidy
ref <- "Sb313" 

FileName <- paste0(paste(spec,collapse="-"), ".tiff") #creates file name for plot image files

#direction structure follows: /rawGenomes/[speciesID]/[versionID]/annotation
#from https://github.com/jtlovell/GENESPACE
#"When working with your own data, place the raw annotation files in this same directory structure with separate directories for each species, separate subdirectories for each genome version, and the annotation files in a subdirectory called "annotation"."

##This is the directory where the GENESPACE will run
#runwd <- file.path("/work/LAS/mhufford-lab/snodgras/Fractionation/GENESPACE_prelim") #Version 0 ran here originally
#for the version 1 update, converting old results to new format 

###MAKE SURE THE RIGHT DIRECTORY IS SPECIFIED!!!###
wd<- "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2"
runwd<-"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2"
genomeRepo<-"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2/genomeRepo"
#Cannot clean the sorghum gene names if this step is to work
#this step makes the directories "bed" and "peptide"
parsed_paths<-parse_annotations(rawGenomeRepo = genomeRepo, 
                                genomeDirs = spec, 
                                genespaceWd = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2",
                                headerEntryIndex = 1,
                                gffIdColumn = "ID")
#see line 11 for how these parameters (spec, ploidy) are set
#Want to make sure that the n unique sequences is equal to the n matched to gff

#this makes the directories: dotplots, pangenes, results, riparian, syntenicHits, tmp
gpar <- init_genespace(
  genomeIDs = spec,
  ploidy = ploidy,
  wd = runwd,
  #diamondUltraSens = TRUE, #(using default for prelim.3)
  orthofinderInBlk = TRUE, 
  nCores = 64,
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = "/opt/MCScanX")

#USING THIS AFTER THE UPDATE TO LATEST GENESPACE
out <- run_genespace(gpar, overwrite = F)


#Customizing output riparian plot
#Customizing the order of the genomes:
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  #c("#5E0D59", "#8A1465", "#B61B5E", "#E12346", "#ED524A", "#EE9E7C","#FFB185","#FC9E40","#FEB820","#FFD53D") #purple, orange, yellow
  #c("#001219","#005F73","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#BB3E03","#AE2012","#9B2226") #Blue to red
  c("#d60000","#f16803","#ECC901","#9DEC01","#00d07f","#03dfc3","#0389df","#031cdf","#7f03df","#CC03DF") #bright rainbow
  )
png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/edited.allvall.riparian.Sb313ref.png", width = 8.5, height = 11, units = "in", res = 600)
p<-plot_riparian(gsParam = out,
                 refGenome = "Sb313",
                 palette = customPal,
                 braidAlpha = .75,
                 chrFill = "lightgrey",
                 addThemes = ggthemes,
                 #chrExpand = 0.0001,
                 genomeIDs = c("Av","Sb313","TdFL","ZnPI615697","ZdMomo","ZdGigi","ZhRIMHU001","ZxTIL25","ZxTIL18","ZvTIL01","ZvTIL11","ZmTzi8",    
                               "ZmNC358","ZmNC350","ZmKi3","ZmKi11","ZmCML69","ZmCML52","ZmCML333","ZmCML322","ZmCML277","ZmCML247",  
                               "ZmCML228","ZmCML103","ZmTx303","ZmMo18W","ZmM37W","ZmHP301","ZmP39","ZmOh7b","ZmOh43","ZmMS71",    
                               "ZmM162W","ZmKy21","ZmIL14H","ZmB97","ZmB73"),
                 #customRefChrOrder = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
                 #scalePlotHeight = 3
                 )

dev.off()

png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2/riparian/edited.allvall.subset.riparian.Sb313ref.png", width = 8.5, height = 11, units = "in", res = 600)
q<-plot_riparian(gsParam = out,
                 refGenome = "Sb313",
                 palette = customPal,
                 braidAlpha = .75,
                 chrFill = "lightgrey",
                 addThemes = ggthemes,
                 #chrExpand = 0.0001,
                 genomeIDs = c("Av","Sb313","TdFL","ZnPI615697","ZdMomo","ZdGigi","ZhRIMHU001","ZxTIL25","ZxTIL18","ZvTIL01","ZvTIL11","ZmB73"),
                 #customRefChrOrder = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
                 #scalePlotHeight = 3
)

dev.off()

#roi <- data.frame(genome = c("ZmB73","ZmB73","ZmB73"),
#                  chr = c("chr1","chr5","chr9"),
#                  color = c("#C8102E","#F1BE48","#F1BE48")) #designate a region of interest
#png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.B73chr159.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
#q<-plot_riparian(gsParam = out,
#                 highlightBed = roi,
#                 refGenome = "ZmB73",
#                 backgroundColor = NULL,
#                 genomeIDs = c("Sb313","Av","TdFL","ZmB73","ZmCML333","ZmNC358","ZmOh43","ZvTIL01","ZxTIL18","ZdGigi"),
#                 customRefChrOrder = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"),
                 #gapProp = 0.02,
                 #scaleBraidGap = 2,
#                 scaleGapSize = 0.5)
#dev.off()
#save(q, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.B73chr159.riparian.Sb313ref.png",type="png")

roi <- data.frame(genome = c("Sb313"),
                  chr = c("Chr01"),
                  color = c("#5E0D59")) #designate a region of interest
png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2/riparian/edited.Sb313chr01.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
r<-plot_riparian(gsParam = out,
                 highlightBed = roi,
                 refGenome = "Sb313",
                 braidAlpha = .75,
                 chrFill = "lightgrey",
                 addThemes = ggthemes,
                 backgroundColor = NULL,
                 genomeIDs = c("TdFL","Sb313","ZmB73"),
                 customRefChrOrder = c("Chr01"),
                 scalePlotHeight = 3)
dev.off()
#save(r, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.Sb313chr10.riparian.Sb313ref.png",type = "png")

roi <- data.frame(genome = c("Sb313","Sb313","Sb313","Sb313","Sb313","Sb313","Sb313","Sb313","Sb313","Sb313"),
                  chr = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
                  #color = c("#5E0D59", "#8A1465", "#B61B5E", "#E12346", "#ED524A", "#EE9E7C","#FFB185","#FC9E40","#FEB820","#FFD53D")
                  c("#d60000","#f16803","#ECC901","#9DEC01","#00d07f","#03dfc3","#0389df","#031cdf","#7f03df","#CC03DF") #bright rainbow
                  ) #designate a region of interest
png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2/riparian/edited.Sb313TdFLZmB73.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
r<-plot_riparian(gsParam = out,
                 highlightBed = roi,
                 refGenome = "Sb313",
                 braidAlpha = .75,
                 chrFill = "lightgrey",
                 addThemes = ggthemes,
                 backgroundColor = NULL,
                 genomeIDs = c("TdFL","Sb313","ZmB73"),
                 scalePlotHeight = 3)
dev.off()


#roi <- data.frame(genome = c("Sb313"),
#                  chr = c("Chr1"),
#                  color = c("#C8102E")) #designate a region of interest
#png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.Sb313chr10.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
#r<-plot_riparian(gsParam = out,
#                 highlightBed = roi,
#                 refGenome = "Sb313",
#                 backgroundColor = NULL,
#                 genomeIDs = c("Sb313","Av","TdFL","ZmB73","ZmCML333","ZmNC358","ZmOh43","ZvTIL01","ZxTIL18","ZdGigi"),
#                 customRefChrOrder = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
#                 scalePlotHeight = 3)
#dev.off()
phased_blk<-read_delim(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/Sb313_phasedBlks.csv", 
                       delim = ",",
                       col_names = TRUE) #reads in the phased_block file

png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/edited.4Subset.Chr10.riparian.Sb313ref.png", width = 4.25, height = 5.5, units = "in", res = 300)
#roi <- data.frame(genome = c("ZmB73","ZmB73","ZmB73"),
#                  chr = c("chr5","chr9","chr6"),
#                  color = c("#C8102E","#C8102E","#F1BE48"))#designate a region of interest
roi <- data.frame(genome = c("Sb313"),
                  chr = c("Chr10"),
                  color = c("#F1BE48"))
Zealution.1<-plot_riparian(gsParam = out,
                 highlightBed = roi,
                 refGenome = "Sb313",
                 backgroundColor = NULL,
                 genomeIDs = c("Sb313","TdFL","ZdGigi","ZvTIL01","ZmB73"),
                 #customRefChrOrder = c("Chr01","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"),
                 #gapProp = 0.02,
                 #scaleBraidGap = 2,
                 scaleGapSize = 0.5)
dev.off()

png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/edited.TdFLvsZmB73.riparian.Sb313ref.png", width = 4.25, height = 5.5, units = "in", res = 300)
Zealution.2<-plot_riparian(gsParam = out,
                           refGenome = "ZmB73",
                           backgroundColor = NULL,
                           genomeIDs = c("TdFL","ZmB73")#,
                           #customRefChrOrder = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
                           #scalePlotHeight = 3
                           )
dev.off

filter(phased_blk, refChr == "Chr10" & grepl('ZmB73', blkID)) %>% arrange(.$blkID) %>% View() #see that the blkIDs are duplicated from the reciprocal comparisons of genomes (focusing only on Sb313 chr10 and blocks involving ZmB73)

uniqblkID<-filter(phased_blk, refChr == "Chr10" & grepl('ZmB73', blkID)) %>% select(blkID) %>% unique() #create a vector with their unique BlockID strings

Sb313chr10.ZmB73blks<-filter(phased_blk, refChr == "Chr10" & grepl('ZmB73', blkID)) %>% arrange(.$blkID) #make it into a dataframe
#get rid of the repeated blocks
filtered.Sb313chr10.ZmB73blks<-Sb313chr10.ZmB73blks[match(uniqblkID$blkID,Sb313chr10.ZmB73blks$blkID),] 
filter(filtered.Sb313chr10.ZmB73blks, genome1=="Sb313" | genome2 =="Sb313") %>% View() #THEN LOOK TO SEE IF THE BLOCKS THAT ARE ZMB73 and SB313 CREATE 2 SETS OF B73 COORDINATES FOR THE SAME SPOT ON CHR10 IN SB313 AND WHAT DO YOU KNOW?? THEY DO! 

####################################################################
#Now that we have filtered things down to what we want to look at, how do we separate
library(GenomicRanges)
gr<-filter(filtered.Sb313chr10.ZmB73blks, genome1=="Sb313" | genome2 =="Sb313") %>%
  GRanges(seqnames = Rle(.$chr1),
          ranges = IRanges(start = .$startBp1, end = .$endBp1, names = .$blkID),
          strand = Rle(strand(.$orient)),
          queryChr = .$chr2,
          queryStart = .$startBp2,
          queryEnd = .$endBp2
  )
gr #not sure why it added all other metadata columns, but it works!
gr[,c("blkID","chr2")] #only display the blkID and chr2 metadata columns
s.gr<-split(gr, gr$chr2) #splits the genomic range object by the query chr

intersect(s.gr[1],s.gr[2],ignore.strand=TRUE) #this will only return the x values and no metadata
subsetByOverlaps(s.gr[1],s.gr[2], ignore.strand = TRUE) %>% unlist() #this returns the x values AND the metadata

#so let's have y be the query chr that has the largest total width of all fragments
#We'd assign the y to file 1 and then if there's something returned from x, it'd go into file 2 (because it overlaps)
t<-ranges(s.gr)[2]
width(t) %>% unlist() %>% str()
sum(width(t) %>% unlist()) %>% paste("index","number")
#OK, so now that we can access the widths from the split list (by query chr)
t<-c()
for(i in 1:length(s.gr)){
  t<-c(t,sum(width(ranges(s.gr)[i]) %>% unlist()))
}
t
#need to find the max length
y<-grep(max(t),t)
#loop it through the intersect
for(i in 1:length(s.gr)){
  if(i != y){
    assign(paste("intersect",i,"by",y, sep = "."),subsetByOverlaps(s.gr[i],s.gr[y], ignore.strand = TRUE) %>% unlist())
  }
  else{
    assign("dup.1.set", unlist(s.gr[y]))
  }
}
intersect.1.by.2
intersect.3.by.2
#So the above are the coordinates that need to go into the 2nd file
if(length(s.gr) > 2){ #first check that this is even needed by checking if there are more than 2 query chromosomes for this region
  for(i in 1:length(s.gr)){ 
    for(j in 1:length(s.gr)){
      if(i != y & j != y & j != i){ #make sure i and j != each other and != the largest chunk y
        if(subsetByOverlaps(s.gr[i],s.gr[j],ignore.strand = TRUE) %>% length() != 0){
         assign(paste("intersect",i,"by",j, sep = "."), subsetByOverlaps(s.gr[i],s.gr[j],ignore.strand = TRUE) %>% unlist()) 
        }
        else{print(paste("There is no overlap between",i, "and",j, "indexes of the split gr object"))}
      }
    }
  }
}else{print("There are not more than 2 query chromosomes to look at the overlap for")}
#this works, but it will repeat (i=1 and j=3 is a repeat of i=3 and j=1)

#HOW TO WRITE THE COORDINATES TO THE DIFFERENT FILES?
#library(rtracklayer)
#export.bed(con = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/test.file2.1.bed", object = intersect.1.by.2)
#only exports the sorghum coordinates...

dup.regions.2<-data.frame(seqnames=c(seqnames(intersect.1.by.2),seqnames(intersect.3.by.2)),
                         starts=c(start(intersect.1.by.2)-1,start(intersect.3.by.2)-1), #-1 because bed files are 0 indexed and GR is 1 indexed
                         ends=c(end(intersect.1.by.2),end(intersect.3.by.2)),
                         names=c(intersect.1.by.2$blkID,intersect.3.by.2$blkID),
                         scores=c(rep(".", length(intersect.1.by.2)),rep(".", length(intersect.3.by.2))),
                         strands=c(strand(intersect.1.by.2),strand(intersect.3.by.2)),
                         query_chr=c(intersect.1.by.2$queryChr,intersect.3.by.2$queryChr),
                         query_start=c(intersect.1.by.2$queryStart-1,intersect.3.by.2$queryStart-1),
                         query_end=c(intersect.1.by.2$queryEnd,intersect.3.by.2$queryEnd))

#HOW TO GET THE ONES THAT DIDN'T DUPLICATE
dup.regions.1<-data.frame(seqnames = c(seqnames(dup.1.set)),
                          starts=c(start(dup.1.set)-1),
                          ends=end(dup.1.set),
                          names=mcols(dup.1.set)$blkID,
                          scores=c(rep(".",length(dup.1.set))),
                          strands=strand(dup.1.set),
                          query_chr=mcols(dup.1.set)$queryChr,
                          query_start=mcols(dup.1.set)$queryStart-1,
                          query_end=mcols(dup.1.set)$queryEnd)

#QC checks:
#1. Total n rows of the dup.regions dataframes == nrow of the filtered df
nrow(filter(filtered.Sb313chr10.ZmB73blks, genome1=="Sb313" | genome2 =="Sb313")) == sum(c(nrow(dup.regions.1),nrow(dup.regions.2)))
#2. All blockIDs found in the filtered df are also found in the duplicated regions:
filter(filtered.Sb313chr10.ZmB73blks, genome1=="Sb313" | genome2 =="Sb313") %>% select("blkID") %>%
  filter(!blkID %in% c(dup.regions.1$names,dup.regions.2$names)) %>% nrow() == 0
#basically filters it for any block ID that is NOT in the duplicated regions; if it = 0, then all filtered block IDs are found in the duplicated regions
#3. All blockIDs found in the duplicated region files are in there exactly once
duplicated(dup.regions.1$names) %>% sum() #should sum to 0 if all are FALSE (means no duplicates)
duplicated(dup.regions.2$names) %>% sum()
duplicated(c(dup.regions.1$names,dup.regions.2$names)) %>% sum()


####NOW THAT WE'VE WORKED OUT HOW TO DO IT FOR ONE REF CHROMOSOME, LET'S BUILD A LOOP OVER ALL THE REF CHROMOSOMES FOR 1 QUERY GENOME
#FIRST MAKE A FUNCTION FOR EACH STEP, MAKE SURE THE FUNCTION WORKS
#Step 0: Filter by ref chromosome, query genome name, and query x ref blocks; make into gr object
step0<-function(REFchr,queryGENOME,refGENOME){ #REFchr = "Chr10", queryGENOME = "ZmB73",refGENOME = "Sb313"
  uniqblkID<-filter(phased_blk, refChr == REFchr & grepl(queryGENOME, blkID)) %>% select(blkID) %>% unique() #create a vector with their unique BlockID strings
  df<-filter(phased_blk, refChr == REFchr & grepl(queryGENOME, blkID)) %>% arrange(.$blkID) #make it into a dataframe
  df<-df[match(uniqblkID$blkID,df$blkID),]
  df<-filter(df, genome1 == refGENOME | genome2 == refGENOME)
  return(df)
}
#test with chr10
t1<-step0("Chr10","ZmB73","Sb313")

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
#test with chr10
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
step3<-function(s.df, y){
  for(i in 1:length(s.df)){
    if(i != y){
      assign(paste("int",i,"by",y, sep = "."),subsetByOverlaps(s.df[i],s.df[y], ignore.strand = TRUE) %>% unlist(), envir = parent.frame())
    }
    else{
      assign("dup.1", unlist(s.df[y]), envir = parent.frame())
    }
  }
  if(length(s.df) > 2){ #first check that this is even needed by checking if there are more than 2 query chromosomes for this region
    for(i in 1:length(s.df)){ 
      for(j in 1:length(s.df)){
        if(i != y & j != y & j != i){ #make sure i and j != each other and != the largest chunk y;
          if(subsetByOverlaps(s.df[i],s.df[j],ignore.strand = TRUE) %>% length() != 0){ #if there's an overlap
              assign(paste("int",i,"by",j, sep = "."), subsetByOverlaps(s.df[i],s.df[j],ignore.strand = TRUE) %>% unlist(),envir = parent.frame()) 
            }
          }
          else{return(print(paste("There is no overlap between",i, "and",j, "indexes of the split gr object")))}
      }
    }
  }else{return(print("There are not more than 2 query chromosomes to look at the overlap for"))}  
}
step3(s.t1,y)

#split step 3 into two parts
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

int.1.df<-step4.a(int.1.by.2)
int.3.df<-step4.a(int.3.by.2)
dup.1.df<-step4.a(dup.1)

dup.2<-step4.b(int.1.df,int.3.df)

#Step 5: Quality control the output to make sure there's no errors

step5<-function(orig.df,dup1,dup2){
  QC1<-nrow(orig.df) == sum(c(nrow(dup1),nrow(dup2)))
  QC2<-filter(orig.df) %>% select("blkID") %>%
    filter(!blkID %in% c(dup1$names,dup2$names)) %>% nrow() == 0
  QC3<-c(duplicated(dup1$names) %>% sum(), #should sum to 0 if all are FALSE (means no duplicates)
         duplicated(dup2$names) %>% sum(),
         duplicated(c(dup1$names,dup2$names)) %>% sum()) %>% paste()
  message<-cat(paste("There are the same number of rows between original df & duplicate dfs:",QC1,"\n",
                 "All block IDs found in the orignal df are found in the duplicate dfs:",QC2,"\n",
                 "All block IDs found in dup1, dup2, and dup1&2 exactly once (0,0,0):", QC3))
  return(message)
}
step5(t1,dup.1.df,dup.2)
#let's break this function into 3 parts for easier looping through; original function is good for manual inspection
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
step5.a(t1,dup.1.df,dup.2) #TRUE TO PASS
step5.b(t1,dup.1.df,dup.2) #TRUE TO PASS
step5.c(t1,dup.1.df,dup.2) #0 TO PASS

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
Sb313_vs_ZmB73.dupregions.2<-step6("Sb313_vs_ZmB73.dupregions.2",dup.2)

#SECOND: Loop through the 10 chromosomes:
remove(list = c(ls(pattern = ".by."),"dup.1","dup.1.df","dup.2.df","t1","s.t1","y"))
for(i in c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09")){#,"Chr10")){
  for(j in c("ZmB73")){ #for when we're eventually looping through the different genomes
    t1<-step0(i,j,"Sb313") #Step 0: Filter by ref chromosome, query genome name, and query x ref blocks; make into gr object
    s.t1<-step1(t1) #Step 1: Split gr object by query chromosomes
    y<-step2(s.t1) #Step 2: Find the query chromosome that contains the most ref bps (max width)
    step3.a(s.t1,y) #Step 3: Overlap the ranges of the smaller width chromosomes with the widest chr from 2
    step3.b(s.t1,y)
    #Step 4: Save the gr objects from 3 as temp dataframes with specific column names to easily make a bed file
    for(k in ls(pattern=".by.")){
      assign(paste0(k,".df",sep=""), step4.a(get(k))) #convert to df each intersect file
    }
    if(length(ls(pattern = ".by.[0-9].df")) > 1){ #if there's multiple intersect dfs, combine them
      for(m in ls(pattern = ".by.[0-9].df")){
      if(exists("dup.2.df")){
        dup.2.df<-step4.b(dup.2.df,get(m)) #if the dataframe already exists, add to it
      }
        else{dup.2.df<-get(m)} #if it doesn't exist, make it with the first int.df
      }
    }
    else{dup.2.df <-get(ls(pattern = ".by.[0-9].df"))} #if there's not multiple intersect dfs
    dup.1.df<-step4.a(dup.1)
    #Step 5: QC to pass
    if(step5.a(t1,dup.1.df,dup.2.df) & #TRUE TO PASS
       step5.b(t1,dup.1.df,dup.2.df) &#TRUE TO PASS
       step5.c(t1,dup.1.df,dup.2.df) == 0 ){#0 TO PASS
      #Step 6: writes the dup.regions to the final files 
      assign(paste0("Sb313_vs_",j,"dupregions.1"), step6(paste0("Sb313_vs_",j,".dupregions.1"), dup.1.df)) 
      assign(paste0("Sb313_vs_",j,"dupregions.2"), step6(paste0("Sb313_vs_",j,".dupregions.2"), dup.2.df))
    }
    else{print(paste("Sb313",i,"by",j,"did not pass QC; NOT WRITTEN TO FINAL FILES"))}
    remove(list = c(ls(pattern = ".by."),"dup.1","dup.1.df","dup.2.df","t1","s.t1","y")) #remove the objects with the same names before the next round of the loop (for genomes)
  }
}
#But doesn't solve the issue of there being an unresolved triplicated region for Chr05
#There are also more than 2 small chunks for Chr07 and Chr08, so will need to figure out how to deal with that for Step3.b
#Step 6 is overwriting the previous information... not adding to it...

t1<-step0("Chr07","ZmB73","Sb313")
s.t1<-step1(t1)
y<-step2(s.t1)
step3.a(s.t1,y)
step3.b(s.t1,y)
for(k in ls(pattern=".by.")){
  assign(paste0(k,".df",sep=""), step4.a(get(k))) #convert to df each intersect file
}
dup1.df<-step4.a(dup.1)

t2<-step0("Chr08","ZmB73","Sb313")
t3<-step0("Chr05","ZmB73","Sb313")
t4<-step0("Chr10","ZmB73","Sb313")

#####
#Let's try finding the coverage/read depth of the Sorghum coordinates with the query genome mapping to it
#To use `coverage()`
t1.cvg<-coverage(s.t1) %>% unlist()
str(t1.cvg)
#t1.cvg<-subset(t1.cvg, t1.cvg <= 2 & t1.cvg > 0) #filters out the values that are 3+ and where there's no coverage; but I don't think it's quite right
gr<-t1 %>%
  GRanges(seqnames = Rle(.$chr1),
          ranges = IRanges(start = .$startBp1, end = .$endBp1, names = .$blkID),
          strand = Rle(strand(.$orient)),
          queryChr = .$chr2,
          queryStart = .$startBp2,
          queryEnd = .$endBp2)
if(t1.cvg[1:62301]@values == 1){print("It worked on values!")}
t1.cvg[]@values
length(t1.cvg@values)
s<-c(1) #vector for starts
e<-c() #vector for ends
v<-c() #vector for coverage values
for(i in 1:length(t1.cvg@values)){
  e<-c(e,sum(t1.cvg@lengths[1:i])) #end will be the cumulative lengths
  v<-c(v,t1.cvg@values[i])
}
for(i in 1:length(e)){
  s<-c(s,e[i-1]+1)
}

gr.cvg<-GRanges(seqnames = c(rep("Chr07",40)),
                ranges = IRanges(start = s, end = e),
                coverage = v)
coverage.filter<-c(1,2)
gr.cvg<-gr.cvg[which(elementMetadata(gr.cvg)[,1] %in% coverage.filter)] #filters the gr object for the coverage data
df.cvg<-data.frame(start = s, end = e, coverage = v)
df.cvg<-filter(df.cvg, coverage %in% coverage.filter)

subsetByOverlaps(gr, gr.cvg)
length(gr)

#restrict(gr, start = ranges(gr.cvg)@start, end = ranges(gr.cvg)@start+ranges(gr.cvg)@width)
restrict(gr, start = ranges(gr.cvg)@start, end = ranges(gr.cvg)@start+ranges(gr.cvg)@width) #didn't work
#invalid start length
restrict(gr, start=c(3,1000000), end=c(7,3000000))
restrict(gr, start = pull(df.cvg$start), end = pull(df.cvg$end)) #didn't work
#invalid start length

###THE LAST ENTRY FOR COVERAGE IS OUTSIDE OF THE START BOUNDS OF THE GR
###FIGURE OUT HOW THAT HAPPENED AND THEN WE MIGHT GET THIS COVERAGE THING SOLVED!
t1.cvg@lengths %>% sum()
max(t1$endBp1)
max(t1$startBp1)

if(max(df.cvg$start) > max(t1$startBp1)){print("coverage start is greater than result starts")}
max(df.cvg$start[df.cvg$start != max(df.cvg$start)]) > max(t1$startBp1) #check if the 2nd to largest startbp

restrict(gr, start = df.cvg$start[1:35], end = df.cvg$end[1:35]) #didn't work
restrict(gr, 
         start = c(16130),#78431,275141,285279,2792526,6345194,6400635,9696801,10251745,12657966,14199575,22820416,37052625),
         end = c(78430)#,275140,285278,2728752,6345193,6400634,9696800,10251744,12657965,14199574,22820415,28115479,49783070)
         )
test.restrict<-c()
for(i in 1:nrow(df.cvg)){test.restrict<-append(test.restrict, restrict(gr, start = df.cvg$start[i], end = df.cvg$end[i]))}
test.restrict
names(test.restrict) #lots of repeated blkIDs

split.test.restrict<-split(test.restrict, test.restrict$chr2)
split.test.restrict[[1]] %>% reduce()
split.test.restrict<-GenomicRanges::reduce(split.test.restrict)
split.test.restrict[4] #THIS PUTS US RIGHT BACK WHERE WE STARTED :((((((((
###BUT what if we start with splitting by query chr, reduce, THEN do the coverage, then restrict
s.t1<-GenomicRanges::reduce(s.t1)
coverage(s.t1) %>% unlist() == t1.cvg #even with reduction, we get the same coverage...
restrict(unlist(s.t1), start = df.cvg$start[1], end = df.cvg$end[1])
reduced.gr.cvg<-GenomicRanges::reduce(gr.cvg) 
s<-c(reduced.gr.cvg@ranges@start) #vector for starts
e<-c(reduced.gr.cvg@ranges@start+reduced.gr.cvg@ranges@width-1) #vector for ends

restrict(unlist(s.t1), start = s, end = e)
s.t1$chr6
reduced.gr.cvg
restrict(s.t1$chr6, start = s[1], end = e[1])
restrict(s.t1$chr10, start = s[1:2], end = e[1:2])

###DO we even need to check the Sb313 coordinates for minor overlaps??? OR just remove the noise? 
#Let's practice with Chr08 instead
#get the coverage of Chr08
s.t2<-step1(t2)
s.t2<-GenomicRanges::reduce(s.t2)
t2.cvg<-coverage(s.t2) %>% unlist()
#get the coordinates for those coverage values
s<-c(1) #vector for starts
e<-c() #vector for ends
v<-c() #vector for coverage values
for(i in 1:length(t2.cvg@values)){
  e<-c(e,sum(t2.cvg@lengths[1:i])) #end will be the cumulative lengths
  v<-c(v,t2.cvg@values[i])
}
for(i in 1:length(e)){
  s<-c(s,e[i-1]+1)
}
#filter the coverage values/Sb313 coordinates
gr2.cvg<-GRanges(seqnames = c(rep("Chr08",32)),
                ranges = IRanges(start = s, end = e),
                coverage = v)
gr2.cvg<-gr2.cvg[which(elementMetadata(gr2.cvg)[,1] %in% coverage.filter)] #filters the gr object for the coverage data
#reduce the filtered coverage Sb313 coordinates
reduced.gr2.cvg<-GenomicRanges::reduce(gr2.cvg) 
s<-c(reduced.gr2.cvg@ranges@start) #vector for starts
e<-c(reduced.gr2.cvg@ranges@start+reduced.gr2.cvg@ranges@width-1) #vector for ends
#Let's try an restrict the original coordinates by the filtered reduced coverage Sb313 coordinates
restrict(unlist(s.t2), start = s, end = e)
restrict(s.t2[1], start = s[4], end = e[4])
restrict(s.t2[3], start = s[1], end = e[1])#undoes what the coverage was supposed to remove
#Maybe try without the reduce step for the coverage coordinates? 
s<-c(gr2.cvg@ranges@start)
e<-c(gr2.cvg@ranges@start+gr2.cvg@ranges@width-1)
restrict(s.t2[5], start = s[2], end = e[2])
#max(gr2.cvg@ranges@start[gr2.cvg@ranges@start != max(gr2.cvg@ranges@start)]) > max(t2$startBp1)
restrict(unlist(s.t2), start = s[1:length(s)-1], end = e[1:length(e)-1]) #doesn't work because length(s) > length(s.t2)
s.t2<-step1(t2) #try it without the reduce t2 step; still same length...

####OK... Let's see what the overlap/coordinate calls are like for all the Sb313 chromosomes 
test.01<-step0("Chr01","ZmB73","Sb313")
test.02<-step0("Chr02","ZmB73","Sb313")
test.03<-step0("Chr03","ZmB73","Sb313")
test.04<-step0("Chr04","ZmB73","Sb313")
test.05<-step0("Chr05","ZmB73","Sb313")
test.06<-step0("Chr06","ZmB73","Sb313")
test.07<-step0("Chr07","ZmB73","Sb313")
test.08<-step0("Chr08","ZmB73","Sb313")
test.09<-step0("Chr09","ZmB73","Sb313")
test.10<-step0("Chr10","ZmB73","Sb313")

split.test.01<-step1(test.01)
split.test.02<-step1(test.02)
split.test.03<-step1(test.03)
split.test.04<-step1(test.04)
split.test.05<-step1(test.05)
split.test.06<-step1(test.06)
split.test.07<-step1(test.07)
split.test.08<-step1(test.08)
split.test.09<-step1(test.09)
split.test.10<-step1(test.10)

ZmB73_phased_blk<-read_delim(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/riparian/ZmB73_phasedBlks.csv", 
                       delim = ",",
                       col_names = TRUE) #reads in the phased_block file
step0.a<-function(REFchr,queryGENOME,refGENOME,phaseddf){ #REFchr = "Chr10", queryGENOME = "ZmB73",refGENOME = "Sb313"
  uniqblkID<-filter(phaseddf, refChr == REFchr & grepl(queryGENOME, blkID)) %>% select(blkID) %>% unique() #create a vector with their unique BlockID strings
  df<-filter(phaseddf, refChr == REFchr & grepl(queryGENOME, blkID)) %>% arrange(.$blkID) #make it into a dataframe
  df<-df[match(uniqblkID$blkID,df$blkID),]
  df<-filter(df, genome1 == refGENOME | genome2 == refGENOME)
  return(df)
}
ZmB73.1<-step0.a("chr1","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.1, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.2<-step0.a("chr2","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.2, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.3<-step0.a("chr3","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.3, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.4<-step0.a("chr4","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.4, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.5<-step0.a("chr5","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.5, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.6<-step0.a("chr6","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.6, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.7<-step0.a("chr7","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.7, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.8<-step0.a("chr8","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.8, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.9<-step0.a("chr9","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.9, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()
ZmB73.10<-step0.a("chr10","Sb313","ZmB73",ZmB73_phased_blk)
select(ZmB73.10, c("chr1","chr2","startBp2","endBp2","orient")) %>% View()

#Check if the overlaps of ZmB73 coordinates of the same Sb313 chromosome correspond to neighboring Sb313 coordinates or distant
select(ZmB73.4, c("chr1","startBp1","endBp1","chr2","startBp2","endBp2","orient","blkID")) %>% View()

#Let's double check how this will work with TdFL
Tdtest.01<-step0("Chr01","TdFL","Sb313")
Tdtest.02<-step0("Chr02","TdFL","Sb313")
Tdtest.03<-step0("Chr03","TdFL","Sb313")
Tdtest.04<-step0("Chr04","TdFL","Sb313")
Tdtest.05<-step0("Chr05","TdFL","Sb313")
Tdtest.06<-step0("Chr06","TdFL","Sb313")
Tdtest.07<-step0("Chr07","TdFL","Sb313")
Tdtest.08<-step0("Chr08","TdFL","Sb313")
Tdtest.09<-step0("Chr09","TdFL","Sb313")
Tdtest.10<-step0("Chr10","TdFL","Sb313")
select(Tdtest.10, c("chr1","startBp1","endBp1","chr2","startBp2","endBp2","orient","blkID")) %>% View()

#make a practice set
test.10 #Sb313 vs. ZmB73
Tdtest.10 #Sb313 vs. TdFL
#ZmB73 vs TdFL for only Sb313 Chr10
#Check to see if order matters (it doesn't)
#step0("Chr10","TdFL","ZmB73") == step0("Chr10","ZmB73","TdFL")
Td.B73.test.10<-step0("Chr10","TdFL","ZmB73")
View(select(Td.B73.test.10, c("genome1","chr1","startBp1","endBp1","genome2","chr2","startBp2","endBp2","orient","refChr","refGenome")))
View(select(Tdtest.10, c("genome1","chr1","startBp1","endBp1","genome2","chr2","startBp2","endBp2","orient","refChr","refGenome")))

B73.endfile.1<-select(test.10, c("chr2","startBp2","endBp2","chr1","startBp1","endBp1","orient")) %>% 
  mutate(namecol = paste(chr1,":",startBp1,"-",endBp1, sep = "")) %>% 
  select(-c("chr1","startBp1","endBp1")) %>%
  filter(chr2 %in% c("chr9","chr5"))
B73.endfile.2<-select(test.10, c("chr2","startBp2","endBp2","chr1","startBp1","endBp1","orient")) %>% 
  mutate(namecol = paste(chr1,":",startBp1,"-",endBp1, sep = "")) %>% 
  select(-c("chr1","startBp1","endBp1")) %>%
  filter(chr2 %in% c("chr6"))
TdFL.endfile.1<-select(Tdtest.10, c("chr2","startBp2","endBp2","chr1","startBp1","endBp1","orient")) %>% 
  mutate(namecol = paste(chr1,":",startBp1,"-",endBp1, sep = "")) %>% 
  select(-c("chr1","startBp1","endBp1")) %>%
  filter(chr2 %in% c("chr4","chr10"))
TdFl.endfile.2<-select(Tdtest.10, c("chr2","startBp2","endBp2","chr1","startBp1","endBp1","orient")) %>% 
  mutate(namecol = paste(chr1,":",startBp1,"-",endBp1, sep = "")) %>% 
  select(-c("chr1","startBp1","endBp1")) %>%
  filter(chr2 %in% c("chr13"))
write_tsv(test.10,"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/Sb313_10vsB73.tsv")
write_tsv(Tdtest.10,"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/Sb313_10vsTdFL.tsv")
write_tsv(Td.B73.test.10, "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/B73vsTdFL.onlySb313Chr10.tsv")
write_tsv(B73.endfile.1,"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/B73.endfile.1.tsv")
write_tsv(B73.endfile.2,"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/B73.endfile.2.tsv")
write_tsv(TdFL.endfile.1,"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/TdFL.endfile.1.tsv")
write_tsv(TdFl.endfile.2,"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/test.parsing.files/TdFL.endfile.2.tsv")

