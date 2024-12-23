#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

cd fullrun.1/genomeRepo
#for i in  Td-FL Zd-Gigi Zd-Momo Zh-RIMHU001 Zn-PI615697 Zv-TIL01 Zv-TIL11 Zx-TIL18 Zx-TIL25 Av; do
#	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3 ] 
#	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3.gz 
#	fi
#	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ]
#	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta.gz
#	fi	
#	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ${i//\-}/${i//\-};
#done

#for i in Zm-B73 Zm-B97 Zm-CML103 Zm-CML228 Zm-CML247 Zm-CML277 Zm-CML322 Zm-CML333 Zm-CML52 Zm-CML69 Zm-HP301 Zm-Ki11 Zm-Ki3 Zm-Ky21 Zm-M162W Zm-M37W Zm-Mo18W Zm-NC350 Zm-NC358 Zm-Oh43 Zm-P39 Zm-Tx303 Zm-Tzi8 ; do
#	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3 ] 
#	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3.gz
#	fi
#	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ]
 #       then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta.gz
  #      fi
#	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ${i//\-}/${i//\-} ;
#done


bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-Il14H*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/Zm-Il14H*.fasta ZmIL14H/ZmIL14H 

bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-Ms71*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/Zm-Ms71*.fasta ZmMS71/ZmMS71

bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-Oh7b*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/Zm-Oh7b*.fasta ZmOh7B/ZmOh7B

#for Sorghum
#bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor*.fa Sb313/Sb313
