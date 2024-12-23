#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=740G   # maximum memory per node
#SBATCH --job-name="genespace"

ml singularity

#ml r-tidyverse

singularity exec --bind /work/LAS/mhufford-lab/snodgras/Fractionation genespace_1.3.1.sif Rscript /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/GenespaceInUnix.R
