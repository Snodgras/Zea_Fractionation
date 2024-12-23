#!/bin/bash
#Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=96:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="Sb-trf"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive

ml trf

trf /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa 2 7 7 80 10 50 2000 -l 1 -d -h


