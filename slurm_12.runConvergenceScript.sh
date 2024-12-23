#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=96:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="convergence"
#SBATCH --exclusive # with this you can exclude asking memory (assuming you request entire node) 


ml r-tidyverse

Rscript /work/LAS/mhufford-lab/snodgras/Fractionation/convergence_2024-12/base_convergence_script.R
