#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=02:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --job-name="GVCFparsing"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

mkdir swap.gvcf_parsing
echo "Script 05.A started at..."
date +%F_%T 
for i in *swap.gvcf; do bash scripts/05.A.cleanDelsOnly.sh $i ; mv cleanDELs* swap.gvcf_parsing/. ; done
echo "Script 05.A finished at..."
date +%F_%T 
echo "Script 05.B started at..." 
date +%F_%T 
for i in swap.gvcf_parsing/cleanDELs*.gvcf; do bash scripts/05.B.GVCF2BED.sh $i ; done
echo "Script 05.B finished at..."
date +%F_%T 
