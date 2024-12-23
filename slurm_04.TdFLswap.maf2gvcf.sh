#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="MAF2GVCF"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/elena/.conda/phg_latest.sif /tassel-5-standalone/run_pipeline.pl \
        -Xmx300g \
        -MAFToGVCFPlugin \
        -referenceFasta /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta \
        -mafFile Sb313_TdFL_anchorwave.swap.maf \
        -sampleName Sb313_TdFL \
        -gvcfOutput Sb313_TdFL_anchorwave.swap.gvcf \
        -fillGaps true > TdFL.out_gvcf.log 2>&1