#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml orthofinder
orthofinder -t $SLURM_JOB_CPUS_PER_NODE -M msa -A mafft -T raxml -b /work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/orthofinder/Results_Aug21/WorkingDirectory/
#orthofinder -o /work/LAS/mhufford-lab/snodgras/Fractionation/betterspeciestree/ -t $SLURM_JOB_CPUS_PER_NODE -M msa -A mafft -T raxml -b /work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/orthofinder/Results_Aug21/WorkingDirectory/
#orthofinder -f /work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/peptide -o /work/LAS/mhufford-lab/snodgras/Fractionation/betterspeciestree/ -t $SLURM_JOB_CPUS_PER_NODE -M msa -A mafft -T raxml -b /work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/orthofinder/Results_Aug21/WorkingDirectory/
