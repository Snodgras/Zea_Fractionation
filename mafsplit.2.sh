#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=05:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="MAF2GVCF"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

ml python/2.7.18-2ut3ogj #or however python2 is invoked in Nova, if necessary
#ml kentutils #or whichever module has this tool in Nova
ml singularity

for input in swap_*.maf
  do
        echo $input
        output=$(echo ${input} | sed 's/.maf//g; s/swap_//g')
        echo $output

cat ${input} | python /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/anchorwave-maf-swap.py  > ${output}_split.maf


#STEP 4:
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName ${output}/ ${output}_split.maf

done

for f in */*.maf ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp"_"$fl";

done

