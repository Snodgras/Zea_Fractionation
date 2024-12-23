#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --time=6:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml samtools
cd /work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_maf2bam/

echo Genome Sb_cnt_aligned_sites Sb_cnt_matched_sites Exon_cnt_aligned_sites Exon_cnt_matched_sites > depth.metrics.txt
for i in *.bam ; do
	sbd=$(samtools depth $i | wc -l ) #sorghum depth of aligned regions
	sbm=$(samtools depth $i | awk '$3>0 {print $0}' | wc -l ) #sorghum depth of matched regions
	ed=$(awk -v OFS='\t' '{print "Chr"$0}' ../ref_Sb313.cds.bed | samtools depth $i -b - | wc -l ) #ref exon depth of aligned regions
	em=$(awk -v OFS='\t' '{print "Chr"$0}' ../ref_Sb313.cds.bed | samtools depth $i -b - | awk '$3>0{print $0}' | wc -l ) #ref exon depth of matched regions
	g=$(echo $i | cut -f 2 -d "_")
	echo $g $sbd $sbm $ed $em >> depth.metrics.txt
done
tr " " "\t" depth.metrics.txt > depth.metrics.tsv
