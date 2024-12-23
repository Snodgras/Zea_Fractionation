#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=03:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load bedtools2

bash scripts/splitGenomes.sh ZdMomo.bin1.fasta assemblies_final/Zd-Momo-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZdMomo.bin1.bed
bash scripts/splitGenomes.sh ZhRIMHU001.bin1.fasta assemblies_final/Zh-RIMHU001-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZhRIMHU001.bin1.bed
bash scripts/splitGenomes.sh ZnPI615697.bin1.fasta assemblies_final/Zn-PI615697-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZnPI615697.bin1.bed
bash scripts/splitGenomes.sh ZvTIL01.bin1.fasta assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZvTIL01.bin1.bed
bash scripts/splitGenomes.sh ZvTIL11.bin1.fasta assemblies_final/Zv-TIL11-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZvTIL11.bin1.bed
bash scripts/splitGenomes.sh ZxTIL18.bin1.fasta assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZxTIL18.bin1.bed
bash scripts/splitGenomes.sh ZxTIL25.bin1.fasta assemblies_final/Zx-TIL25-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZxTIL25.bin1.bed
bash scripts/splitGenomes.sh ZmB73.bin1.fasta NAM-assemblies/Zm-B73.fasta parsing_coordinates/zmB73vsZmB73.bin1.bed
bash scripts/splitGenomes.sh ZmB97.bin1.fasta NAM-assemblies/Zm-B97.fasta parsing_coordinates/zmB73vsZmB97.bin1.bed
bash scripts/splitGenomes.sh ZmCML103.bin1.fasta NAM-assemblies/Zm-CML103.fasta parsing_coordinates/zmB73vsZmCML103.bin1.bed
bash scripts/splitGenomes.sh ZmCML228.bin1.fasta NAM-assemblies/Zm-CML228.fasta parsing_coordinates/zmB73vsZmCML228.bin1.bed
bash scripts/splitGenomes.sh ZmCML247.bin1.fasta NAM-assemblies/Zm-CML247.fasta parsing_coordinates/zmB73vsZmCML247.bin1.bed
bash scripts/splitGenomes.sh ZmCML277.bin1.fasta NAM-assemblies/Zm-CML277.fasta parsing_coordinates/zmB73vsZmCML277.bin1.bed
bash scripts/splitGenomes.sh ZmCML322.bin1.fasta NAM-assemblies/Zm-CML322.fasta parsing_coordinates/zmB73vsZmCML322.bin1.bed
bash scripts/splitGenomes.sh ZmCML333.bin1.fasta NAM-assemblies/Zm-CML333.fasta parsing_coordinates/zmB73vsZmCML333.bin1.bed
bash scripts/splitGenomes.sh ZmCML52.bin1.fasta NAM-assemblies/Zm-CML52.fasta parsing_coordinates/zmB73vsZmCML52.bin1.bed
bash scripts/splitGenomes.sh ZmCML69.bin1.fasta NAM-assemblies/Zm-CML69.fasta parsing_coordinates/zmB73vsZmCML69.bin1.bed
bash scripts/splitGenomes.sh ZmHP301.bin1.fasta NAM-assemblies/Zm-HP301.fasta parsing_coordinates/zmB73vsZmHP301.bin1.bed
bash scripts/splitGenomes.sh ZmIL14H.bin1.fasta NAM-assemblies/Zm-IL14H.fasta parsing_coordinates/zmB73vsZmIL14H.bin1.bed
bash scripts/splitGenomes.sh ZmKi11.bin1.fasta NAM-assemblies/Zm-Ki11.fasta parsing_coordinates/zmB73vsZmKi11.bin1.bed
bash scripts/splitGenomes.sh ZmKi3.bin1.fasta NAM-assemblies/Zm-Ki3.fasta parsing_coordinates/zmB73vsZmKi3.bin1.bed
bash scripts/splitGenomes.sh ZmKy21.bin1.fasta NAM-assemblies/Zm-Ky21.fasta parsing_coordinates/zmB73vsZmKy21.bin1.bed
bash scripts/splitGenomes.sh ZmM162W.bin1.fasta NAM-assemblies/Zm-M162W.fasta parsing_coordinates/zmB73vsZmM162W.bin1.bed
bash scripts/splitGenomes.sh ZmM37W.bin1.fasta NAM-assemblies/Zm-M37W.fasta parsing_coordinates/zmB73vsZmM37W.bin1.bed
bash scripts/splitGenomes.sh ZmMo18W.bin1.fasta NAM-assemblies/Zm-Mo18W.fasta parsing_coordinates/zmB73vsZmMo18W.bin1.bed
bash scripts/splitGenomes.sh ZmMS71.bin1.fasta NAM-assemblies/Zm-MS71.fasta parsing_coordinates/zmB73vsZmMS71.bin1.bed
bash scripts/splitGenomes.sh ZmNC350.bin1.fasta NAM-assemblies/Zm-NC350.fasta parsing_coordinates/zmB73vsZmNC350.bin1.bed
bash scripts/splitGenomes.sh ZmNC358.bin1.fasta NAM-assemblies/Zm-NC358.fasta parsing_coordinates/zmB73vsZmNC358.bin1.bed
bash scripts/splitGenomes.sh ZmOh43.bin1.fasta NAM-assemblies/Zm-Oh43.fasta parsing_coordinates/zmB73vsZmOh43.bin1.bed
bash scripts/splitGenomes.sh ZmOh7b.bin1.fasta NAM-assemblies/Zm-Oh7b.fasta parsing_coordinates/zmB73vsZmOh7b.bin1.bed
bash scripts/splitGenomes.sh ZmP39.bin1.fasta NAM-assemblies/Zm-P39.fasta parsing_coordinates/zmB73vsZmP39.bin1.bed
bash scripts/splitGenomes.sh ZmTx303.bin1.fasta NAM-assemblies/Zm-Tx303.fasta parsing_coordinates/zmB73vsZmTx303.bin1.bed
bash scripts/splitGenomes.sh ZmTzi8.bin1.fasta NAM-assemblies/Zm-Tzi8.fasta parsing_coordinates/zmB73vsZmTzi8.bin1.bed

bash scripts/splitGenomes.sh ZdMomo.bin2.fasta assemblies_final/Zd-Momo-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZdMomo.bin2.bed
bash scripts/splitGenomes.sh ZhRIMHU001.bin2.fasta assemblies_final/Zh-RIMHU001-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZhRIMHU001.bin2.bed
bash scripts/splitGenomes.sh ZnPI615697.bin2.fasta assemblies_final/Zn-PI615697-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZnPI615697.bin2.bed
bash scripts/splitGenomes.sh ZvTIL01.bin2.fasta assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZvTIL01.bin2.bed
bash scripts/splitGenomes.sh ZvTIL11.bin2.fasta assemblies_final/Zv-TIL11-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZvTIL11.bin2.bed
bash scripts/splitGenomes.sh ZxTIL18.bin2.fasta assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZxTIL18.bin2.bed
bash scripts/splitGenomes.sh ZxTIL25.bin2.fasta assemblies_final/Zx-TIL25-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/zmB73vsZxTIL25.bin2.bed
bash scripts/splitGenomes.sh ZmB73.bin2.fasta NAM-assemblies/Zm-B73.fasta parsing_coordinates/zmB73vsZmB73.bin2.bed
bash scripts/splitGenomes.sh ZmB97.bin2.fasta NAM-assemblies/Zm-B97.fasta parsing_coordinates/zmB73vsZmB97.bin2.bed
bash scripts/splitGenomes.sh ZmCML103.bin2.fasta NAM-assemblies/Zm-CML103.fasta parsing_coordinates/zmB73vsZmCML103.bin2.bed
bash scripts/splitGenomes.sh ZmCML228.bin2.fasta NAM-assemblies/Zm-CML228.fasta parsing_coordinates/zmB73vsZmCML228.bin2.bed
bash scripts/splitGenomes.sh ZmCML247.bin2.fasta NAM-assemblies/Zm-CML247.fasta parsing_coordinates/zmB73vsZmCML247.bin2.bed
bash scripts/splitGenomes.sh ZmCML277.bin2.fasta NAM-assemblies/Zm-CML277.fasta parsing_coordinates/zmB73vsZmCML277.bin2.bed
bash scripts/splitGenomes.sh ZmCML322.bin2.fasta NAM-assemblies/Zm-CML322.fasta parsing_coordinates/zmB73vsZmCML322.bin2.bed
bash scripts/splitGenomes.sh ZmCML333.bin2.fasta NAM-assemblies/Zm-CML333.fasta parsing_coordinates/zmB73vsZmCML333.bin2.bed
bash scripts/splitGenomes.sh ZmCML52.bin2.fasta NAM-assemblies/Zm-CML52.fasta parsing_coordinates/zmB73vsZmCML52.bin2.bed
bash scripts/splitGenomes.sh ZmCML69.bin2.fasta NAM-assemblies/Zm-CML69.fasta parsing_coordinates/zmB73vsZmCML69.bin2.bed
bash scripts/splitGenomes.sh ZmHP301.bin2.fasta NAM-assemblies/Zm-HP301.fasta parsing_coordinates/zmB73vsZmHP301.bin2.bed
bash scripts/splitGenomes.sh ZmIL14H.bin2.fasta NAM-assemblies/Zm-IL14H.fasta parsing_coordinates/zmB73vsZmIL14H.bin2.bed
bash scripts/splitGenomes.sh ZmKi11.bin2.fasta NAM-assemblies/Zm-Ki11.fasta parsing_coordinates/zmB73vsZmKi11.bin2.bed
bash scripts/splitGenomes.sh ZmKi3.bin2.fasta NAM-assemblies/Zm-Ki3.fasta parsing_coordinates/zmB73vsZmKi3.bin2.bed
bash scripts/splitGenomes.sh ZmKy21.bin2.fasta NAM-assemblies/Zm-Ky21.fasta parsing_coordinates/zmB73vsZmKy21.bin2.bed
bash scripts/splitGenomes.sh ZmM162W.bin2.fasta NAM-assemblies/Zm-M162W.fasta parsing_coordinates/zmB73vsZmM162W.bin2.bed
bash scripts/splitGenomes.sh ZmM37W.bin2.fasta NAM-assemblies/Zm-M37W.fasta parsing_coordinates/zmB73vsZmM37W.bin2.bed
bash scripts/splitGenomes.sh ZmMo18W.bin2.fasta NAM-assemblies/Zm-Mo18W.fasta parsing_coordinates/zmB73vsZmMo18W.bin2.bed
bash scripts/splitGenomes.sh ZmMS71.bin2.fasta NAM-assemblies/Zm-MS71.fasta parsing_coordinates/zmB73vsZmMS71.bin2.bed
bash scripts/splitGenomes.sh ZmNC350.bin2.fasta NAM-assemblies/Zm-NC350.fasta parsing_coordinates/zmB73vsZmNC350.bin2.bed
bash scripts/splitGenomes.sh ZmNC358.bin2.fasta NAM-assemblies/Zm-NC358.fasta parsing_coordinates/zmB73vsZmNC358.bin2.bed
bash scripts/splitGenomes.sh ZmOh43.bin2.fasta NAM-assemblies/Zm-Oh43.fasta parsing_coordinates/zmB73vsZmOh43.bin2.bed
bash scripts/splitGenomes.sh ZmOh7b.bin2.fasta NAM-assemblies/Zm-Oh7b.fasta parsing_coordinates/zmB73vsZmOh7b.bin2.bed
bash scripts/splitGenomes.sh ZmP39.bin2.fasta NAM-assemblies/Zm-P39.fasta parsing_coordinates/zmB73vsZmP39.bin2.bed
bash scripts/splitGenomes.sh ZmTx303.bin2.fasta NAM-assemblies/Zm-Tx303.fasta parsing_coordinates/zmB73vsZmTx303.bin2.bed
bash scripts/splitGenomes.sh ZmTzi8.bin2.fasta NAM-assemblies/Zm-Tzi8.fasta parsing_coordinates/zmB73vsZmTzi8.bin2.bed
