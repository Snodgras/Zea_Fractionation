#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=48:00:00
#SBATCH --job-name=gatk-chr6
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml gatk
ref=/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.chr6.fasta
ml samtools
# index
#samtools faidx $ref
# dict
#gatk CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${ref%.*}.dict
# db import
gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
-V ready_Sb313_TdFL_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_TdKS_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdGigi_4to1_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdMomo_4to1_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZhRIMHU001_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB73_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB97_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML103_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML228_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML247_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML277_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML322_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML333_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML52_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML69_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmHP301_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmIL14H_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmKi11_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmKi3_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmKy21_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmM162W_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmM37W_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmMo18W_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmMS71_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmNC350_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmNC358_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmOh43_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmOh7b_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmP39_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmTx303_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmTzi8_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZnPI615697_4to1_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL01_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL11_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZxTIL18_Chr06.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZxTIL25_Chr06.bin2.gvcf.gz.recode.vcf \
--batch-size 1 \
--genomicsdb-workspace-path chr6.bin2_gatkDBimport \
-L 06 \
--genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 
#--tmp-dir $TMPDIR

# vcf output
gatk --java-options "-Xmx50g" GenotypeGVCFs \
-R $ref \
-V gendb://chr6.bin2_gatkDBimport \
--cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000 \
-O chr6.bin2.vcf