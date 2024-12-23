#libraries
library(tidyverse)

#Read in files
ID.key<-read.table("Sb.coord.ID.key", header = T,stringsAsFactors = F)

#The "#" at the beginning of the header line will cause issues for reading in the file 
#The lack of "END" will be an issue too
Zl_RIL003.gvcf<-read.table("delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.test",col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Zl-RIL003-Reference-PanAnd-2.0","ASM_Chr","ASM_End","ASM_Start","ASM_Strand","Sb.coord.ID"),stringsAsFactors = F)
Zv_TIL01.gvcf<-read.table("delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.test",col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Zv-TIL01-Reference-PanAnd-2.0","ASM_Chr","ASM_End","ASM_Start","ASM_Strand","Sb.coord.ID"),stringsAsFactors = F)
Zv_TIL11.gvcf<-read.table("delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.test",col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Zv-TIL11-Reference-PanAnd-2.0","ASM_Chr","ASM_End","ASM_Start","ASM_Strand","Sb.coord.ID"),stringsAsFactors = F)
Zx_TIL18.gvcf<-read.table("delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.test",col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Zx-TIL18-Reference-PanAnd-2.0","ASM_Chr","ASM_End","ASM_Start","ASM_Strand","Sb.coord.ID"),stringsAsFactors = F)
Zx_TIL25.gvcf<-read.table("delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.test",col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Zx-TIL25-Reference-PanAnd-2.0","ASM_Chr","ASM_End","ASM_Start","ASM_Strand","Sb.coord.ID"),stringsAsFactors = F)

#Reformat files
Zv_TIL11.gvcf<-mutate(Zv_TIL11.gvcf, GT = str_split(string = Zv_TIL11.gvcf[,10], pattern = ":", simplify = TRUE)[,1])
new.Zv_TIL11<-select(Zv_TIL11.gvcf, c("Sb.coord.ID","ASM_Chr","ASM_Start","ASM_End","ASM_Strand","CHROM","POS","ALT","REF","GT")) %>% arrange(.,Sb.coord.ID)
duplicated(new.Zv_TIL11$Sb.coord.ID) %>% summary()
Zv_TIL11.dupIDs<-filter(new.Zv_TIL11, duplicated(new.Zv_TIL11$Sb.coord.ID) == TRUE) %>% select("Sb.coord.ID") %>% pull() 
filter(new.Zv_TIL11, Sb.coord.ID %in% Zv_TIL11.dupIDs)

collapsed.Zv_TIL11<-tibble(Sb.coord.ID=NA,
                           ASM_Chr=NA,
                           ASM_Start=NA,
                           ASM_End=NA,
                           ASM_Strand=NA,
                           CHROM=NA,
                           POS=NA,
                           ALT=NA,
                           REF=NA,
                           GT=NA)
for(i in 1:nrow(new.Zv_TIL11)){
  if(new.Zv_TIL11$Sb.coord.ID[i] %in% Zv_TIL11.dupIDs & !new.Zv_TIL11$Sb.coord.ID[i] %in% collapsed.Zv_TIL11$Sb.coord.ID){
    tmp<-filter(new.Zv_TIL11, Sb.coord.ID == new.Zv_TIL11$Sb.coord.ID[i])
    collapsed.Zv_TIL11<-add_row(collapsed.Zv_TIL11,
                                Sb.coord.ID=new.Zv_TIL11$Sb.coord.ID[i],
                                ASM_Chr=toString(tmp$ASM_CHR),
                                ASM_Start=toString(tmp$ASM_Start),
                                ASM_End=toString(tmp$ASM_End),
                                ASM_Strand=toString(tmp$ASM_Strand),
                                CHROM=toString(tmp$CHROM),
                                POS=toString(tmp$POS),
                                ALT=toString(tmp$ALT),
                                REF=toString(tmp$REF),
                                GT=toString(tmp$GT))
  } else{
    if(!new.Zv_TIL11$Sb.coord.ID[i] %in% collapsed.Zv_TIL11$Sb.coord.ID){
      collapsed.Zv_TIL11<-add_row(collapsed.Zv_TIL11,
                                  Sb.coord.ID=new.Zv_TIL11$Sb.coord.ID[i],
                                  ASM_Chr=new.Zv_TIL11$ASM_Chr[i],
                                  ASM_Start=new.Zv_TIL11$ASM_Start[i],
                                  ASM_End=new.Zv_TIL11$ASM_End[i],
                                  ASM_Strand=new.Zv_TIL11$ASM_Strand[i],
                                  CHROM=new.Zv_TIL11$CHROM[i],
                                  POS=new.Zv_TIL11$POS[i],
                                  ALT=new.Zv_TIL11$ALT[i],
                                  REF=new.Zv_TIL11$REF[i],
                                  GT=new.Zv_TIL11$GT[i])
    }
  }
}
collapsed.Zv_TIL11<-na.omit(collapsed.Zv_TIL11)

reformat.gvcf<-function(df){
  df<-mutate(df, GT = str_split(string = df[,10], pattern = ":", simplify = TRUE)[,1]) #makes genotype a column
  df<-select(df, c("Sb.coord.ID","ASM_Chr","ASM_Start","ASM_End","ASM_Strand","CHROM","POS","ALT","REF","GT")) %>% arrange(.,Sb.coord.ID) #selects columns and sorts by Sb ID
  df.dupIDs<-filter(df, duplicated(df$Sb.coord.ID) == TRUE) %>% select("Sb.coord.ID") %>% pull() #creates a list of duplicate IDs
  collapsed.df<-tibble(Sb.coord.ID=NA,
                             ASM_Chr=NA,
                             ASM_Start=NA,
                             ASM_End=NA,
                             ASM_Strand=NA,
                             CHROM=NA,
                             POS=NA,
                             ALT=NA,
                             REF=NA,
                             GT=NA) #initiates an empty df for values to be added to
  for(i in 1:nrow(df)){ 
    #if the SB id is in the duplicate list and hasn't already been added to the collapsed df
    if(df$Sb.coord.ID[i] %in% df.dupIDs & !df$Sb.coord.ID[i] %in% collapsed.df$Sb.coord.ID){
      tmp<-filter(df, Sb.coord.ID == df$Sb.coord.ID[i]) #create df with only values of duplicated id
      #Add row, columns of duplicate df (tmp) are collapsed into a string (value, value) by "toString"
      collapsed.df<-add_row(collapsed.df,
                                  Sb.coord.ID=df$Sb.coord.ID[i],
                                  ASM_Chr=toString(tmp$ASM_CHR),
                                  ASM_Start=toString(tmp$ASM_Start),
                                  ASM_End=toString(tmp$ASM_End),
                                  ASM_Strand=toString(tmp$ASM_Strand),
                                  CHROM=toString(tmp$CHROM),
                                  POS=toString(tmp$POS),
                                  ALT=toString(tmp$ALT),
                                  REF=toString(tmp$REF),
                                  GT=toString(tmp$GT))
    } else{ #Otherwise if the ID hasn't been added to the collapsed df yet
      if(!df$Sb.coord.ID[i] %in% collapsed.df$Sb.coord.ID){
        #add it as a row to the collapsed df (no collapse needed since it isn't a duplicate)
        collapsed.df<-add_row(collapsed.df,
                                    Sb.coord.ID=df$Sb.coord.ID[i],
                                    ASM_Chr=df$ASM_Chr[i],
                                    ASM_Start=df$ASM_Start[i],
                                    ASM_End=df$ASM_End[i],
                                    ASM_Strand=df$ASM_Strand[i],
                                    CHROM=df$CHROM[i],
                                    POS=df$POS[i],
                                    ALT=df$ALT[i],
                                    REF=df$REF[i],
                                    GT=df$GT[i])
      }
    }
  }
  collapsed.df<-na.omit(collapsed.df) #remove initializing NA row
  return(collapsed.df) #Return the collapsed df
}
collapsed.Zl_RIL003<-reformat.gvcf(Zl_RIL003.gvcf)
collapsed.Zv_TIL01<-reformat.gvcf(Zv_TIL01.gvcf)
collapsed.Zx_TIL18<-reformat.gvcf(Zx_TIL18.gvcf)
collapsed.Zx_TIL25<-reformat.gvcf(Zx_TIL25.gvcf)
