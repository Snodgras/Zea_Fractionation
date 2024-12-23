#!/bin/bash
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --mem=350GB
#SBATCH -t 4:00:00
#SBATCH -J featureCounts.
#SBATCH -o featureCounts.o%j
#SBATCH -e featureCounts.e%j
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

for g in ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	while read -r query sb bin ; do 
		if [ -f "$g"_"$sb"."$bin".maf ] ; then 
			tail -n +2 "$g"_"$query"_"$sb"_sorted_filtered.maf >> "$g"_"$sb"."$bin".maf
			else
			cat "$g"_"$query"_"$sb"_sorted_filtered.maf > "$g"_"$sb"."$bin".maf
		fi
	done < /work/LAS/mhufford-lab/snodgras/Fractionation/parsing_coordinates/${g}.bins.txt
done

for g in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
for i in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 ; do
	mv ${g}_${i}.bin1*.maf ${g}_${i}.bin1.maf
	mv ${g}_${i}.bin2*.maf ${g}_${i}.bin2.maf
done
done
