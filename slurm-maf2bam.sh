#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --time=6:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml samtools
cd /work/LAS/mhufford-lab/snodgras/Fractionation/
for i in AnchorWave_output/Sb313_Z*anchorwave.maf ; do 
	python maf-convert.py sam ${i} | samtools view -O BAM --reference Sb313.fasta - | samtools sort - >  AnchorWave_maf2bam/${i#AnchorWave_output/}.bam
done

