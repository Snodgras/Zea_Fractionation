#!/usr/bin/env Rscript
#Running Genespace in unix instead of R Studio
#ml singularity
#SINGULARITY_IMAGE=genespace_1.3.1.sif

library(GENESPACE)

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
namesdf <-  names[names$spec %in% reqNames ,]

spec <- namesdf$spec
ploidy <- namesdf$ploidy
ref <- "Sb313" 
wd<- "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2"
runwd<-"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2"
genomeRepo<-"/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2/genomeRepo"
parsed_paths<-parse_annotations(rawGenomeRepo = genomeRepo, 
                                genomeDirs = spec, 
                                genespaceWd = "/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.2",
                                headerEntryIndex = 1,
                                gffIdColumn = "ID")
gpar <- init_genespace(
  genomeIDs = spec,
  ploidy = ploidy,
  wd = runwd,
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
roi <- data.frame(genome = c("Sb313"),
                  chr = c("Chr01"),
                  color = c("#d60000")) #designate a region of interest
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
