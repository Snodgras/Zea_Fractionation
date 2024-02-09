library(GENESPACE)
#spec <- c("zTIL01", "zTIL11", "zTIL18", "zTIL25", "zB73v5", "zdgigi", "zdmomo", "zmhuet", "avirgi", "ppanic", "sbicol", "rtuber", "tdacts", "irugos" )
#ploidy <- c(2,2,2,2,2,2,2,2,1,1,1,1,2,1)
#ref <- "sbicol"

#names is a tsv file with the first field the species name (abbrev.) and second field is ploidy
#WILL WE NEED TO SPECIFY THE SORGHUM AS DIPLOID AND EVERYTHING ELSE AS TETRAPLOID (YES)
names <- read.delim("/rawGenomes/names.tsv") #change as needed
reqNames <- c("zTIL11", "zB73v5")
#Specifies the genomes we want to plot so that it speeds up plotting
#Can specify anything we want
#To run full pairs, then we'd use the full table in the namesdf
namesdf <-  names[names$spec %in% reqNames ,]

spec <- namesdf$spec
ploidy <- namesdf$ploidy
ref <- "zB73v5" #WILL I NEED TO CHANGE THIS REF NAME

FileName <- paste0(paste(spec,collapse="-"), ".tiff") #creates file name for plot image files

#direction structure follows: /rawGenomes/[speciesID]/[versionID]/annotation
#from https://github.com/jtlovell/GENESPACE
#"When working with your own data, place the raw annotation files in this same directory structure with separate directories for each species, separate subdirectories for each genome version, and the annotation files in a subdirectory called "annotation"."

##This is the directory where the GENESPACE will run
runwd <- file.path("/work/LAS/mhufford-lab/snodgras/Fractionation/GENESPACE_prelim")
#see line 11 for how these parameters (spec, ploidy) are set
#DOUBLE CHECK ORTHOFINDER SETTINGS
gpar <- init_genespace(
  genomeIDs = spec,
  speciesIDs = spec,
  versionIDs = spec,
  outgroup = NULL,
  ploidy = ploidy,
  diamondMode = "fast",
  orthofinderMethod = "fast",
  wd = runwd,
  orthofinderInBlk = TRUE, 
  overwrite = F, 
  verbose = T,
  nCores = 64,
  minPepLen = 50,
  gffString = "gff",
  pepString = "pep",
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = "/opt/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))
#WILL `pars_annotations()` throw errors because of differences between PanAnd, NAM, and sorghum genomes?
#LOOK AT ARUN'S GFF FILES TO SEE HOW THE FILES ARE FORMATTED (Probably only sorghum needs changing)
parse_annotations(
  gsParam = gpar,
  gffEntryType = "gene",
  gffIdColumn = "ID",
  gffStripText = "ID=",
  headerEntryIndex = 1,
  headerSep = " ",
  headerStripText = "ID=", troubleshoot = TRUE)
#Runs orthofinder
gpar <- run_orthofinder(
  gsParam = gpar)
#Runs synteny (MCScanX??)
gpar <- synteny(gsParam = gpar)
#creates the plot as a tiff file
tiff(
  FileName,
  units = "in",
  width = 13,
  height = 8,
  res = 600
)
#creates the riparian plots
plot_riparianHits(
  gpar,
  chrLabCex = 1,
  genomeLabCex = 1,
  refGenome = ref,
  gapProp = .001,
  blackBg = FALSE,
  verbose = T
)
dev.off()

#saves plot as a tiff file
tiff(
  paste0("ogs_", FileName),
  units = "in",
  width = 13,
  height = 8,
  res = 600
)
#old version of the riparian plot command
plot_riparian(
  gpar,
  chrLabCex = 1,
  genomeLabCex = 1,
  refGenome = ref,
  gapProp = .001,
  blackBg = FALSE,
  verbose = T
)
dev.off()
#WILL WE NEED TO BUILD THE PANGENOME ANNOTATIONS `pg()`???
#SEE STEP 9: https://htmlpreview.github.io/?https://github.com/jtlovell/GENESPACE/blob/master/doc/genespaceOverview.html
#We need to find the pairs of syntenic blocks and the (sorghum) genes that it contains
#Then we need to find which pairs of blocks have duplicate sorghum genes
#Block ID | coordinates on all the different genomes | array of sorghum genes in that block