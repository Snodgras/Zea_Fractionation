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

ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave gff2seq \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \
	-o Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa

