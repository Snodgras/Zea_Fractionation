#!/bin/bash
#Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=1:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --partition=nova

#for i in assemblies_final/*.fasta ; do 
#	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/02.minimapAnchorPoints.sh ${i} 
#done
#for i in NAM-assemblies/*.fasta ; do 
#	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/02.minimapAnchorPoints.sh ${i}
#done

bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/02.minimapAnchorPoints.sh assemblies_final/Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fasta
