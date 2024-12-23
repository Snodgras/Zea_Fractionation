#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=8:00:00

for i in ready_Sb313_*Chr01*bin1*vcf ; do
	echo $i >> problem_indels.txt
	awk -v OFS='\t' 'length($4) > 9101263 || length($5) - 10 > 9101263 {print $2, "ref length:"length($4), "alt length:"length($5)-10}' ${i} >> problem_indels.txt
done

#for j in 02 03 04 05 06 07 08 09 10 ; 
#	do for i in ready_Sb313_*Chr${j}*bin1*.vcf ; 
#		do echo $i >> problem_indels.txt ; 
#		awk -v OFS='\t' 'length($4) > 9101263 || length($5) - 10 > 9101263 {print $2, "ref length:"length($4), "alt length:"length($5)-10}' ${i} >> problem_indels.txt ; 
#	done ; 
#done  

#for j in 01 02 03 04 05 06 07 08 09 10 ; 
#	do for i in ready_Sb313_*Chr${j}*bin2*.vcf ; 
#		do echo $i >> problem_indels.txt ; 
#		awk -v OFS='\t' 'length($4) > 9101263 || length($5) - 10 > 9101263 {print $2, "ref length:"length($4), "alt length:"length($5)-10}' ${i} >> problem_indels.txt ; 
#	done ; 
#done  

#for i in ready_Sb313_*vcf ; do 
#	awk -v OFS='\t' 'length($4) < 9101263 || length($5) - 10 < 9101263 {print $0}' ${i} > pass_${i}
#done
