#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=4:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
ml gatk
#gatk SelectVariants -V chr10.vcf -O chr10.indelonly.vcf --select-type-to-include INDEL

gatk SelectVariants -V chr10.bin1.vcf -O chr10.bin1.indelonly.vcf --select-type-to-include INDEL
gatk SelectVariants -V chr10.bin2.vcf -O chr10.bin2.indelonly.vcf --select-type-to-include INDEL

ml bcftools
bcftools annotate -x INFO,^FORMAT/GT chr10.bin1.indelonly.vcf | grep "^##" -v - > chr10.bin1.headerless.indelonly.vcf
bcftools annotate -x INFO,^FORMAT/GT chr10.bin2.indelonly.vcf | grep "^##" -v - > chr10.bin2.headerless.indelonly.vcf

