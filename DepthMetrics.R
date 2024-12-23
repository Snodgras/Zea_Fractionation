#Depth metrics
library(tidyverse)

depth_metrics<-read_tsv(file = "depth.metrics.tsv")

#Number of bases in Sorghum genome: 732152042
#Number of bases in Ref Exons: 16390866

depth_metrics<-depth_metrics %>% mutate(Sb_percent_aligned = Sb_cnt_aligned_sites/732152042,
                         Sb_percent_matched = Sb_cnt_matched_sites/732152042,
                         Exon_percent_aligned = Exon_cnt_aligned_sites/16390866,
                         Exon_percent_matched = Exon_cnt_matched_sites/16390866)
select(depth_metrics, contains("percent"))

ggplot(depth_metrics, aes(x = Sb_percent_aligned, y = Exon_percent_aligned))+
  geom_point()
ggsave("depth.genomeVsExonAligned.percentage.pdf", device="pdf",dpi=300)

filter(depth_metrics, Genome != "TdFL") %>% 
  ggplot(aes(x = Sb_percent_aligned, y = Exon_percent_aligned))+
  geom_point()

filter(depth_metrics, Genome != "TdFL") %>% 
  ggplot(aes(x = Sb_percent_aligned, y = Exon_percent_aligned))+
  geom_point(alpha = 0.5, color="green")+
  geom_text(aes(label = Genome),size=12/.pt)
