#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=1:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#ml gatk
#gatk VariantsToTable -V chr6.indelonly.vcf -O chr6.indelonly.table -F CHROM -F POS -F TYPE -F NCALLED -GF GT -GF GQ

ml bcftools
bcftools annotate -x INFO,^FORMAT/GT chr6.bin1.indelonly.vcf  > chr6.bin1.indelonly_reformatted.vcf
bcftools annotate -x INFO,^FORMAT/GT chr6.bin2.indelonly.vcf  > chr6.bin2.indelonly_reformatted.vcf
