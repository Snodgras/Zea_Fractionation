library(GENESPACE)
#spec <- c("zTIL01", "zTIL11", "zTIL18", "zTIL25", "zB73v5", "zdgigi", "zdmomo", "zmhuet", "avirgi", "ppanic", "sbicol", "rtuber", "tdacts", "irugos" )
#ploidy <- c(2,2,2,2,2,2,2,2,1,1,1,1,2,1)
#ref <- "sbicol"

#names is a tsv file with the first field the species name (abbrev.) and second field is ploidy
#WILL WE NEED TO SPECIFY THE SORGHUM AS DIPLOID AND EVERYTHING ELSE AS TETRAPLOID (YES)
names <- read.delim("/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/names.tsv")
reqNames <- c("Sb313","TdFL","ZdGigi","ZvTIL01", "ZxTIL18","ZmB73","ZmCML333","ZmNC358","ZmOh43","Av")
#Specifies the genomes we want to plot so that it speeds up plotting
#Can specify anything we want
#To run full pairs, then we'd use the full table in the namesdf
namesdf <-  names[names$spec %in% reqNames ,]

spec <- namesdf$spec
ploidy <- namesdf$ploidy
ref <- "Sb313" #WILL I NEED TO CHANGE THIS REF NAME

FileName <- paste0(paste(spec,collapse="-"), ".tiff") #creates file name for plot image files

#direction structure follows: /rawGenomes/[speciesID]/[versionID]/annotation
#from https://github.com/jtlovell/GENESPACE
#"When working with your own data, place the raw annotation files in this same directory structure with separate directories for each species, separate subdirectories for each genome version, and the annotation files in a subdirectory called "annotation"."

##This is the directory where the GENESPACE will run
#runwd <- file.path("/work/LAS/mhufford-lab/snodgras/Fractionation/GENESPACE_prelim") #Version 0 ran here originally
#for the version 1 update, converting old results to new format 
wd<- "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3"
runwd<-"/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3"
genomeRepo<-"/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/genomeRepo"
#Cannot clean the sorghum gene names if this step is to work
parsed_paths<-parse_annotations(rawGenomeRepo = genomeRepo, 
                                genomeDirs = spec, 
                                genespaceWd = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3",
                                headerEntryIndex = 1,
                                gffIdColumn = "ID")
#see line 11 for how these parameters (spec, ploidy) are set

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
out <- run_genespace(gpar)

#Customizing output riparian plot
#Customizing the order of the genomes:
png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.allvall.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
p<-plot_riparian(gsParam = out,
              refGenome = "Sb313",
              backgroundColor = NULL,
              genomeIDs = c("Sb313","Av","TdFL","ZmB73","ZmCML333","ZmNC358","ZmOh43","ZvTIL01","ZxTIL18","ZdGigi"),
              customRefChrOrder = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
              scalePlotHeight = 3)
dev.off()
#save(p, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.allvall.riparian.Sb313ref.png", type = "png")

roi <- data.frame(genome = c("ZmB73","ZmB73","ZmB73"),
                  chr = c("chr1","chr5","chr9"),
                  color = c("#C8102E","#F1BE48","#F1BE48")) #designate a region of interest
png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.B73chr159.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
q<-plot_riparian(gsParam = out,
              highlightBed = roi,
              refGenome = "ZmB73",
              backgroundColor = NULL,
              genomeIDs = c("Sb313","Av","TdFL","ZmB73","ZmCML333","ZmNC358","ZmOh43","ZvTIL01","ZxTIL18","ZdGigi"),
              customRefChrOrder = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10"),
              #gapProp = 0.02,
              #scaleBraidGap = 2,
              scaleGapSize = 0.5)
dev.off()
#save(q, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.B73chr159.riparian.Sb313ref.png",type="png")

roi <- data.frame(genome = c("Sb313"),
                  chr = c("Chr10"),
                  color = c("#C8102E")) #designate a region of interest
png(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.Sb313chr10.riparian.Sb313ref.png", width = 7, height = 7, units = "in", res = 300)
r<-plot_riparian(gsParam = out,
              highlightBed = roi,
              refGenome = "Sb313",
              backgroundColor = NULL,
              genomeIDs = c("Sb313","Av","TdFL","ZmB73","ZmCML333","ZmNC358","ZmOh43","ZvTIL01","ZxTIL18","ZdGigi"),
              customRefChrOrder = c("Chr01","Chr02","Chr03","Chr04","Chr05","Chr06","Chr07","Chr08","Chr09","Chr10"),
              scalePlotHeight = 3)
dev.off()
#save(r, file = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/edited.Sb313chr10.riparian.Sb313ref.png",type = "png")

#Try seeing how many blocks are pulled out for a given roi
#trying with chromosome 10 of sorghum (because it's the smallest)
qreturn<-query_pangenes(gsParam = out, bed = roi)

library(tidyverse)
SynBlockCoords<-read_csv(file="/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/results/syntenicBlock_coordinates.csv",col_names =TRUE)
str(SynBlockCoords)
unique(SynBlockCoords$blkID)
colnames(SynBlockCoords)

query_pangenes(gsParam = out, bed = tibble(genome="ZmB73",chr="chr9",start=0,end=1000000))

#Want to pull out each syntenic block from the others and make sure it's 1:2 sorghum to everything and 1:1 everything else
#BUT HOW DO I DO THIS? 

#What is the different types of blocks (selfblk vs. other)? How are blocks named? 
#pull out 1 sorghum block and figure out what other pairwise blocks overlap it? Seems like a nightmare...

####From Parsing Genespace Results R script I wrote on my local computer:
##Import phased block csv
phased_blk<-read_delim(file = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/Sb313_phasedBlks.csv", 
                       delim = ",",
                       col_names = TRUE)
##Split block ID variable into front half and back half
phased_blk<-mutate(phased_blk, numOnly_blkID = str_split(blkID, ":", simplify = TRUE)[,2],
                   genOnly_blkID = str_split(blkID, ":", simplify = TRUE)[,1],
                   blklength1 = endBp1-startBp1,
                   blklength2 = endBp2-startBp2) 

summary(select(phased_blk, c(blklength1, blklength2))) #not sure why block 2 has start bps ahead of end bps...?

max(phased_blk$numOnly_blkID) #there are only 99 numbered block ids
### if a block is "4_1" for one genome pair, is it the same for a different genome pair? 
filter(phased_blk, numOnly_blkID == " 4_1") %>% select(c(contains("genome"), contains("chr"), 
                                                         contains("blk"), contains("start"), contains("end"))) %>% View()

#filter(phased_blk, numOnly_blkID == " 221_1") %>% select(c(contains("genome"), contains("chr"), 
#                                                           contains("blk"), contains("start"), contains("end"))) %>% View()

filter(phased_blk, numOnly_blkID == " 66_1") %>% select(c(contains("genome"), contains("chr"), 
                                                          contains("blk"), contains("start"), contains("end"))) %>% View()

#filter(phased_blk, numOnly_blkID == " 325_1") %>% select(c(contains("genome"), contains("chr"), 
#                                                           contains("blk"), contains("start"), contains("end"))) %>% View()

filter(phased_blk, numOnly_blkID == " 1_1") %>% select(c(contains("genome"), contains("chr"), 
                                                         contains("blk"), contains("start"), contains("end"))) %>% View()

filter(phased_blk, numOnly_blkID == " 1_2") %>% select(c(contains("genome"), contains("chr"), 
                                                         contains("blk"), contains("start"), contains("end"))) %>% View()

###Is there a better way to look at this instead of manually?
filter(phased_blk, numOnly_blkID  == " 1_1") %>%
  count(., genOnly_blkID ) %>% 
  #print(n=110)
  right_join(x=filter(phased_blk, numOnly_blkID  == " 1_1"), y=., by="genOnly_blkID") %>%
  ggplot(aes(x = genome1, y=genome2))+
  geom_tile(aes(fill=n))+
  scale_fill_gradient(low="white", high="blue")

filter(phased_blk, numOnly_blkID  == " 1_1") %>%
  count(., genOnly_blkID ) %>% 
  #print(n=110)
  right_join(x=filter(phased_blk, numOnly_blkID  == " 1_1"), y=., by="genOnly_blkID") %>%
  ggplot(aes(x = genome1, y=genome2))+
  geom_tile(aes(fill=refChr))

plot_blkNumber<-function(ID){
  plt<-filter(phased_blk, numOnly_blkID  == ID) %>%
    count(., genOnly_blkID ) %>% 
    right_join(x=filter(phased_blk, numOnly_blkID  == ID), y=., by="genOnly_blkID") %>%
    ggplot(aes(x = genome1, y=genome2))+
    geom_tile(aes(fill=n))+
    scale_fill_gradient(low="yellow", high="blue")+
    ggtitle(paste(ID, " by number of blocks per pair"))
  return(plt)
}

plot_refChr<-function(ID){
  plt<-filter(phased_blk, numOnly_blkID  == ID) %>%
    count(., genOnly_blkID ) %>% 
    right_join(x=filter(phased_blk, numOnly_blkID  == ID), y=., by="genOnly_blkID") %>%
    ggplot(aes(x = genome1, y=genome2))+
    geom_tile(aes(fill=refChr))+
    ggtitle(paste(ID, " by Sb313 Ref Chr per pair"))
  return(plt)
  
}

plot_refChr(" 1_1")
plot_blkNumber(" 1_1")

s<-sample(phased_blk$numOnly_blkID, 3, replace = FALSE)

plot_refChr(s[1])
plot_blkNumber(s[1])

plot_refChr(s[2])
plot_blkNumber(s[2])

plot_refChr(s[3])
plot_blkNumber(s[3])

unique.numOnly_blkID<-unique(phased_blk$numOnly_blkID)
unique.numOnly_blkID[unique.numOnly_blkID == ""]<-"unassigned"

library(here)
for(i in 1:length(unique.numOnly_blkID)){
  x<-str_remove_all(string = unique.numOnly_blkID[i], pattern = " ") #gets rid of extra spaces in the ID name
  plot_refChr(unique.numOnly_blkID[i]) #custom function to plot
  pathName="/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/"
  fileName1=paste0("RefChr.heatmap.",x,".png")
  ggsave(path=pathName, filename = fileName1, device = "png", width = 7, height = 7, units = "in", dpi = 300)
  plot_blkNumber(unique.numOnly_blkID[i]) #second custom function to plot
  fileName2=paste0("BlkNumber.heatmap.",x,".png")
  ggsave(path=pathName, filename = fileName2, device = "png", width = 7, height = 7, units = "in", dpi = 300)
  print(paste("We're ", (i/length(unique.numOnly_blkID))*100, "% complete")) #sanity marker for me to see how far along its running
}

#What if we filter by size of block?
#What about filtering by core/near-core genes? What's the lower threshold for GENESPACE to make a block?
#filter out the orthogroups that are large (gene families?)

#to show the average number of paired blocks for each numeric block ID
phased_blk %>% 
  group_by(numOnly_blkID) %>% 
  count(., genOnly_blkID) %>% summarize(meanPairedBlkCount=mean(n))%>%
  ggplot(aes(x=meanPairedBlkCount))+
  geom_histogram(binwidth = 1)
ggsave(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/avgPairedBlksperBlockID.histogram.png",device="png")

#to show the avergae number of refChromosomes for each numeric block ID
phased_blk %>%
  group_by(numOnly_blkID, refChr) %>% count() %>%
  ggplot(aes(x=refChr, y=numOnly_blkID))+
  geom_tile(aes(fill=n))+
  scale_fill_gradient(low="yellow", high="blue")+
  theme(axis.text.y = element_blank())
ggsave(filename = "/work/LAS/mhufford-lab/snodgras/Fractionation/prelim.3/riparian/countOfRefChrPerBlockID.heatmap.png",device="png")
