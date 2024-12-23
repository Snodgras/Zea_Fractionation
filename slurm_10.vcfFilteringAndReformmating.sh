#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=1:00:00

for i in {2..9} ; do
	for m in bin1 bin2 ; do 
		bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/10.vcfFilteringAndReformatting.sh chr${i}.${m}.headerless.indelonly.vcf 0${i}
		echo Done with ${i} ${m}
	done
done
