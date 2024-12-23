#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --time=4:00:00
#SBATCH --job-name=Chr04.maf2fasta
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err

ml singularity
singularity run https://depot.galaxyproject.org/singularity/phast:1.5--hec16e2b_5

while read line ; do
	echo $line | tr " " "\t" > temp.Chr04.bed
	for g in TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
		maf_parse --features temp.Chr04.bed ${g}_Chr04.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Chr04.bin1.${g}.maf
		maf_parse --features temp.Chr04.bed ${g}_Chr04.bin2.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Chr04.bin2.${g}.maf
		if [$g == TdFL] ; then
			awk -v id=$(cut -f 4 temp.Chr04.bed) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.Chr04.bin1.${g}.maf |tr "\t" "\n" >> Chr04.bin1.fasta
			awk -v id=$(cut -f 4 temp.Chr04.bed) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.Chr04.bin2.${g}.maf |tr "\t" "\n" >> Chr04.bin2.fasta
			else
				awk -v id=$(cut -f 4 temp.Chr04.bed) -v OFS='\t' '$1=="s" & $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Chr04.bin1.${g}.maf |tr "\t" "\n" >> Chr04.bin1.fasta
				awk -v id=$(cut -f 4 temp.Chr04.bed) -v OFS='\t' '$1=="s" & $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Chr04.bin2.${g}.maf |tr "\t" "\n" >> Chr04.bin2.fasta
		fi
	done
	rm temp.Chr04*
done < /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.04.bed
