#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=72:00:00
#SBATCH --job-name=gatk-trim
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

#load module
module load samtools
module load gatk

sed -e 's/Chr//g' /ptmp/LAS/snodgras/Fractionation/Sb313.fasta > Sb313.clean.fasta
samtools faidx Sb313.clean.fasta
gatk CreateSequenceDictionary -R Sb313.clean.fasta
