# dN/dS using Yin et al 2022 MBE data

#Subset out the gene IDs that overlap between our 12K reference genes and their Sorghum genes

library(tidyverse)

ref_gene_list<-read_tsv("/work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.geneids.txt", col_names = c("ExonCnt","GeneID_Sb313"))
Yin2022MBE_data<-readxl::read_excel("/work/LAS/mhufford-lab/snodgras/Fractionation/reformatted-yin2022MBE-suppTable3.xlsx")

dn_ds_data<-filter(Yin2022MBE_data, str_detect(GeneID_Sb313, pattern = paste(str_remove_all(ref_gene_list$GeneID_Sb313,".[1-9].v3.1"), collapse = "|")))

nrow(dn_ds_data) #should be ~3,195

####Let's compare dn/ds between M1 and M2 overall####

mutate(dn_ds_data, w_diff = ω_M1 - ω_M2) %>% #neg diff means dn/ds is higher for M2 than M1, pos diff means dn/ds is higher for M1 than M2
	ggplot(aes(x=w_diff))+
	geom_freqpoly(binwidth = 0.01)+
	theme_bw()

ggplot(dn_ds_data)+
  geom_histogram(binwidth = 0.01, aes(x=ω_M1), fill = "dodgerblue3", alpha =0.5)+
  geom_histogram(binwidth = 0.01, aes(x=ω_M2), fill = "darkred", alpha = 0.5)+
  theme_bw()
	
mutate(dn_ds_data, w_diff = ω_M1 - ω_M2) %>% summarize(w.mean=mean(w_diff,na.rm = T),
                                                       w.median=median(w_diff,na.rm=T),
                                                       w.sd=sd(w_diff,na.rm=T))

t.test(x=dn_ds_data$ω_M1,y=dn_ds_data$ω_M2, paired = T, alternative = "less") #choosing lesser because M2 dN/dS should be larger than M1
#t = -1.9772, df = 3193, p-value = 0.02405
#mean difference: -0.009333469 

# are there any instances where M2 dn/ds was < M1 dn/ds
filter(dn_ds_data, ω_M2 < ω_M1) %>% nrow()/nrow(dn_ds_data)
#Yes, about 46.3% of genes

#### Compare to Deletion overlap ####
long_collapsed_dels$size_class<-factor(long_collapsed_dels$size_class, levels = c("1-10bp","11-100bp","101-1000bp","1001-10Kbp","10K-100Kbp","100K+bp"))
mutate(long_collapsed_dels, Gene_ID.ref = str_remove_all(Gene_ID.ref, ".[1-9].v3.1")) %>%
  inner_join(x=., y=dn_ds_data, by=c("Gene_ID.ref"="GeneID_Sb313")) %>%
  ggplot()+
  geom_histogram(binwidth = 0.01, aes(x=ω_M1), fill = "dodgerblue3", alpha =0.5)+
  geom_histogram(binwidth = 0.01, aes(x=ω_M2), fill = "darkred", alpha = 0.5)+
  theme_bw()+
  xlab("w")+
  facet_grid(rows=vars(size_class))


#### Compare to convergence status ####
temp<-convergence.shared %>% filter(Convergence_Category %in% c("CompletelyShared","CompletelyDifferent","SomeShared"))
temp<-temp %>% mutate(Gene_ID = str_remove_all(Gene_ID,".[1-9]*.v3.[1-9]*"))
temp<-select(dn_ds_data, GeneID_Sb313, ω_M1, ω_M2) %>% left_join(x=., y=temp, by = c("GeneID_Sb313" = "Gene_ID"))

#how many genes in the dn/ds data set are in each category
na.omit(temp) %>% select(GeneID_Sb313, Subgenome, Convergence_Category)  %>% unique() %>% group_by(Subgenome, Convergence_Category) %>% count()
#CompletelyDifferent: M1= 97, M2=72
#SomeShared: M1=723, M2=643
#CompletelyShared: M1=2393, M2=1922

#is dn/ds different for genes that are called completely shared, completely different, some shared in pairwise comparisons? 
library(ggridges)
temp$Convergence_Category<-temp$Convergence_Category %>% factor(levels = c("CompletelyDifferent","SomeShared","CompletelyShared"))
filter(temp, Subgenome == "M1") %>%
  ggplot(aes(x = ω_M1, y= Target_Genome))+
  facet_wrap(vars(Convergence_Category))+
  geom_density_ridges(aes(fill=Target_Genome), show.legend = F)+
  scale_fill_manual(values=genome_colors)+
  theme_bw()+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/dn_ds_ConvergenceCategories.M1.ridgeline.png",
       device="png",dpi= 300, width = 8, height=8)  
filter(temp, Subgenome == "M2") %>%
  ggplot(aes(x = ω_M2, y= Target_Genome))+
  facet_wrap(vars(Convergence_Category))+
  geom_density_ridges(aes(fill=Target_Genome), show.legend = F)+
  scale_fill_manual(values=genome_colors)+
  theme_bw()+ylab("")
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/dn_ds_ConvergenceCategories.M2.ridgeline.png",
       device="png",dpi= 300, width = 8, height=8)  

summary(aov(data = filter(temp, Convergence_Category == "CompletelyDifferent" & Subgenome == "M1"), ω_M1 ~ Target_Genome))
TukeyHSD(aov(data = filter(temp, Convergence_Category == "CompletelyDifferent"& Subgenome == "M1"), ω_M1 ~ Target_Genome))$Target_Genome %>%
  as_tibble(rownames = NA) %>% rownames_to_column()%>% 
  filter(`p adj` <=0.05)

summary(aov(data = filter(temp, Convergence_Category == "SomeShared" & Subgenome == "M1"), ω_M1 ~ Target_Genome))
TukeyHSD(aov(data = filter(temp, Convergence_Category == "SomeShared"& Subgenome == "M1"), ω_M1 ~ Target_Genome))$Target_Genome %>%
  as_tibble(rownames = NA) %>% rownames_to_column()%>% 
  filter(`p adj` <=0.05)

summary(aov(data = filter(temp, Convergence_Category == "CompletelyShared" & Subgenome == "M1"), ω_M1 ~ Target_Genome))
TukeyHSD(aov(data = filter(temp, Convergence_Category == "CompletelyShared" & Subgenome == "M1"), ω_M1 ~ Target_Genome))$Target_Genome %>%
  as_tibble(rownames = NA) %>% rownames_to_column()%>% 
  filter(`p adj` <=0.05)

ggplot(filter(temp, Subgenome == "M1"), 
       aes(x=Target_Genome, y=Query_Genome))+
  geom_tile(aes(fill = ω_M1))+
  facet_wrap(vars(Convergence_Category))+
  theme_bw()+theme(axis.text.x = element_text(angle=90))+
  ggtitle("dN/dS for M1 genes")+xlab("")+ylab("")+
  scale_fill_viridis_c()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/dn_ds_ConvergenceCategories.M1.heatmap.png",
       dpi=300, device="png",width=10,height=8)

ggplot(filter(temp, Subgenome == "M2"), 
       aes(x=Target_Genome, y=Query_Genome))+
  geom_tile(aes(fill = ω_M2))+
  facet_wrap(vars(Convergence_Category))+
  theme_bw()+theme(axis.text.x = element_text(angle=90))+
  ggtitle("dN/dS for M2 genes")+xlab("")+ylab("")+
  scale_fill_viridis_c()
ggsave("/work/LAS/mhufford-lab/snodgras/Fractionation/Fractionation_Plots/dn_ds_ConvergenceCategories.M2.heatmap.png",
       dpi=300, device="png",width=10,height=8)

#### Compare to Fractionation Status####
str(long.full.fractionation.status)
