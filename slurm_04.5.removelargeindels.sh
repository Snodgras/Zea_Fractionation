#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=24:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml vcftools
#for i in Sb313*.gvcf.gz ; do 
#	awk '$0 ~ /Indel is too long/  { print }' ${i}.trimGVCF.cmds_*.e* |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
#	vcftools --gzvcf trimmed_${i} --exclude-positions positions.txt --recode --recode-INFO-all --out ready_${i}
#done
ml gatk
#for i in ready*recode.vcf ; do gatk IndexFeatureFile -I $i ;done

for i in GVCF/*TdKS*.gvcf.gz GVCF/*4to1*.gvcf.gz ; do 
	awk '$0 ~ /Indel is too long/  { print }' ${i#GVCF/}.TdKSZnZd.trimGVCF.cmds_*.e* |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
	vcftools --gzvcf trimmed_${i#GVCF/} --exclude-positions positions.txt --recode --recode-INFO-all --out ready_${i#GVCF/}
done

for i in ready*TdKS*recode.vcf ready*4to1*recode.vcf ; do gatk IndexFeatureFile -I $i ; done
