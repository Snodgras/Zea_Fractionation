#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=12:00:00
#SBATCH --job-name=gatk-chr100
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml gatk
ref=/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.chr100.fasta
ml samtools
# index
samtools faidx $ref
# dict
gatk CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${ref%.*}.dict
# db import
gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
-V ready_Sb313_TdFL_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdGigi_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB73_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL10_Chr10.bin1.gvcf.gz.recode.vcf \
--batch-size 1 \
--genomicsdb-workspace-path chr100.bin1_gatkDBimport \
-L 10 \
--genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 
#--tmp-dir $TMPDIR

# vcf output
gatk --java-options "-Xmx50g" GenotypeGVCFs \
-R $ref \
-V gendb://chr100.bin1_gatkDBimport \
--cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000 \
-O chr100.bin1.vcf

# db import
gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
-V ready_Sb313_TdFL_Chr10.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdGigi_Chr10.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB73_Chr10.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL10_Chr10.bin2.gvcf.gz.recode.vcf \
--batch-size 1 \
--genomicsdb-workspace-path chr100.bin2_gatkDBimport \
-L 10 \
--genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 
#--tmp-dir $TMPDIR

# vcf output
gatk --java-options "-Xmx50g" GenotypeGVCFs \
-R $ref \
-V gendb://chr100.bin2_gatkDBimport \
--cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000 \
-O chr100.bin2.vcf
