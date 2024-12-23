####Libraries to load####
library(tidyverse)

####Handy vectors####

tripsacinae_genomes<-c("TdKS","TdFL","ZnPI615697","ZdGigi","ZdMomo","ZnPI615697_4to1","ZdGigi_4to1","ZdMomo_4to1",
						"ZhRIMHU001","ZxTIL18","ZxTIL25","ZvTIL01","ZvTIL11",
						...)
genome_colors

chr_colors

####inputs to load####

#Sorghum CDS references

#Headerless VCFs per Sb ref chr per subgenome

for(i in 1:10){
	for(b in c("bin1", "bin2")){
		assign(paste0("Chr",i,"_",b,"_vcf"), 
		read_tsv(file = paste0("/work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/chr",i,".",b,".headerless.indelonly.vcf"), 
		colnames = T)
		)
	}
}

####Split alt alleles####

#split into different columns
for(i in ls(***list the names of the vcf objects***)){
	assign(i,
		mutate(i, ALT1 = str_split(ALT, ",",simplify=T)[,1],
			   ALT2 = str_split(ALT, ",",simplify=T)[,2],
			   ALT1_length = case_when(length(ALT1) > 0 ~ length(ALT1) - length(REF),
			   						   length(ALT1) == 0 ~ NA), #double check that this is the correct length command
			   ALT2_length = case_when(length(ALT2) > 0 ~ length(ALT2) - length(REF),
			   						   length(ALT2) == 0 ~ NA), #double check that this is the correct length command
			   ALT1_type = case_when(ALT1_length < 0 ~ "del",
			   						 ALT1_length > 0 ~ "ins",
			   						 is.na(ALT1_length) ~ NA),
			   ALT2_type = case_when(ALT2_length < 0 ~ "del",
			   						 ALT2_length > 0 ~ "ins",
			   						 is.na(ALT2_length) ~ NA)
			   )
	)
}

#split into different rows and filter out insertions and NAs
function(obj)
newdf<-tibble(CHROM = NA, POS=NA, REF=NA, ALT=NA,QUAL=NA,Length=NA,
				TdFL=NA,TdKS=NA, ZdGigi=NA,ZdMomo=NA,ZdGigi_4to1=NA,ZdMomo_4to1=NA,
				ZhRIMHU001=NA,ZmB73=NA,ZmB97=NA,ZmCML103=NA,ZmCML228=NA,ZmCML247=NA,ZmCML277=NA,
				ZmCML322=NA,ZmCML333=NA,ZmCML52=NA,ZmCML69=NA,ZmHP301=NA,ZmIL14H=NA,ZmKi11=NA,
				ZmKi3=NA,ZmNC358=NA,ZmNC350=NA,ZmOh43=NA,ZmOh7b=NA,ZmP39=NA,ZmTx303=NA,ZmTzi8=NA,
				ZvTIL01=NA,ZvTIL11=NA,ZxTIL18=NA,ZxTIL25=NA,ZnPI615697=NA,ZnPI615697_4to1=NA)
for(i in 1:nrow(obj)){
	if(ALT1_type == "del"){
	newdf<-add_row(newdef, CHROM=obj[i,1], POS=obj[i,2], REF=obj[i,4], ALT=obj[i,"ALT1"],QUAL=obj[i,6],Length=obj[i,"ALT1_length"],
				###CHANGE COLUMN NUMBERS###
				TdFL=obj[i,1],
				TdKS=obj[i,1], 
				ZdGigi=obj[i,1],
				ZdMomo=obj[i,1],
				ZdGigi_4to1=obj[i,1],
				ZdMomo_4to1=obj[i,1],
				ZhRIMHU001=obj[i,1],
				ZmB73=obj[i,1],
				ZmB97=obj[i,1],
				ZmCML103=obj[i,1],
				ZmCML228=obj[i,1],
				ZmCML247=obj[i,1],
				ZmCML277=obj[i,1],
				ZmCML322=obj[i,1],
				ZmCML333=obj[i,1],
				ZmCML52=obj[i,1],
				ZmCML69=obj[i,1],
				ZmHP301=obj[i,1],
				ZmIL14H=obj[i,1],
				ZmKi11=obj[i,1],
				ZmKi3=obj[i,1],
				ZmNC358=obj[i,1],
				ZmNC350=obj[i,1],
				ZmOh43=obj[i,1],
				ZmOh7b=obj[i,1],
				ZmP39=obj[i,1],
				ZmTx303=obj[i,1],
				ZmTzi8=obj[i,1],
				ZvTIL01=obj[i,1],
				ZvTIL11=obj[i,1],
				ZxTIL18=obj[i,1],
				ZxTIL25=obj[i,1],
				ZnPI615697=obj[i,1],
				ZnPI615697_4to1=obj[i,1])
	}
	if(ALT2_type == "del"){
	newdf<-add_row(newdef, CHROM=obj[i,1], POS=obj[i,2], REF=obj[i,4], ALT=obj[i,"ALT2"],QUAL=obj[i,6],Length=obj[i,"ALT2_length"],
				###CHANGE COLUMN NUMBERS###
				TdFL=obj[i,1],
				TdKS=obj[i,1], 
				ZdGigi=obj[i,1],
				ZdMomo=obj[i,1],
				ZdGigi_4to1=obj[i,1],
				ZdMomo_4to1=obj[i,1],
				ZhRIMHU001=obj[i,1],
				ZmB73=obj[i,1],
				ZmB97=obj[i,1],
				ZmCML103=obj[i,1],
				ZmCML228=obj[i,1],
				ZmCML247=obj[i,1],
				ZmCML277=obj[i,1],
				ZmCML322=obj[i,1],
				ZmCML333=obj[i,1],
				ZmCML52=obj[i,1],
				ZmCML69=obj[i,1],
				ZmHP301=obj[i,1],
				ZmIL14H=obj[i,1],
				ZmKi11=obj[i,1],
				ZmKi3=obj[i,1],
				ZmNC358=obj[i,1],
				ZmNC350=obj[i,1],
				ZmOh43=obj[i,1],
				ZmOh7b=obj[i,1],
				ZmP39=obj[i,1],
				ZmTx303=obj[i,1],
				ZmTzi8=obj[i,1],
				ZvTIL01=obj[i,1],
				ZvTIL11=obj[i,1],
				ZxTIL18=obj[i,1],
				ZxTIL25=obj[i,1],
				ZnPI615697=obj[i,1],
				ZnPI615697_4to1=obj[i,1])
	}
}

#filter out intergenic only indels