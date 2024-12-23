#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN

ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali -i Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \
	-r Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
	-a Zv-TIL01-Reference-PanAnd-2.0_Sbicolor_313_v3.0_cds.sam \
	-as Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-ar Sbicolor_313_v3.0_ref.sam \
	-s panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta.gz \
	-n align1.anchors \
	-R 2 -Q 1 -ns
