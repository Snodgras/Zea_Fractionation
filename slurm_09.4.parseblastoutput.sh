#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=1:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

cd /work/LAS/mhufford-lab/snodgras/Fractionation/blast_output

for i in *txt ; do 
	awk -v OFS='\t' '$11 < 0.01 {print $0}' ${i} > ${i%.txt}.evalue.filtered
done

for i in *evalue.filtered ; do cut -f 1 $i | sort | uniq -c > ${i}.IDcounts ; done


