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

perl scripts/AlignmentToDotplot.pl Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 Zv-TIL01.CHRONLY_Sbicolor_313_v3.0_cds.CHRONLY.sam > Zv-TIL01.Sbicolor.CHRONLY.tab

