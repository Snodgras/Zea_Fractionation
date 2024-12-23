#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=12:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail


ml bedtools2

#you have to make the GVCF into a bed file with a start and end position (either with the END= info for non-variant entries or the length of the ref allele for variant entries). Otherwise the intersect function will only intersect the points and not the actual blocks within the GVCF
for i in GVCF/Sb313*.gvcf.gz; do
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	o=${i#GVCF/}
	zcat $i | grep -v "#" - | tr ';' '\t' | awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - |sed 's/END=//g' - | sed 's/ASM_Strand=//g' - | bedtools intersect -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | cut -f 4 | sort | uniq > ${o%.gvcf.gz}.aligned.IDs.txt
	zcat $i | grep -v "#" - | tr ';' '\t' | awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - |sed 's/END=//g' - | sed 's/ASM_Strand=//g' - | bedtools intersect -v -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | cut -f 4 | sort | uniq > ${o%.gvcf.gz}.noalignment.IDs.txt
done

#allows each deletion to intersect with the ref exon coordinates
#see if a deletion overlaps with multiple ref exons
#if a del allele is in gvcf, that means the genotype = 1
for i in GVCF/Sb313*.gvcf.gz ; do 
o=${i#GVCF/} #Get the GVCF name
zcat $i | grep -v "#" - | tr ';' '\t' | awk -v OFS="\t" '{print $1,$2,$4,$5,$13}' | sed "s/<NON_REF>/N/g" | awk -v OFS="\t" '{split($4,a,/,/); for(i = 1; i <= length(a); ++i) if(a[i] != "N" && length(a[i]) < length($3)) print $1,$2,$3,$4,a[i],length($3)-length(a[i]),$5}' | awk -v OFS="\t" '{print $1,$2,$2+$6,$3,$4,$5,$6,$7}' | tr ":" "\t" | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${o%.gvcf.gz}.dels.bed
done

for i in *.dels.bed ; do
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	id=$(echo $i | sed 's/.dels.bed/.aligned.IDs.txt/g')
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -a - -b $i > ${i%.dels.bed}.refExons.deleted
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -v -a - -b $i > ${i%.dels.bed}.refExons.retained

	cut -f 4 ${i%.dels.bed}.refExons.deleted | sort | uniq > ${i%.dels.bed}.refExons.deleted.IDs.txt
	cut -f 4 ${i%.dels.bed}.refExons.retained | sort | uniq > ${i%.dels.bed}.refExons.retained.IDs.txt
done

####MAKE SURE THIS IS DONE IN A CLEAN DIRECTORY SO THERE'S NO OVERLAP
for i in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
 echo ${i}.bin1.noAlignment_CDS_ID > ${i}.bin1.allchr.refExons.noalign
 echo ${i}.bin2.noAlignment_CDS_ID > ${i}.bin2.allchr.refExons.noalign
 echo ${i}.bin1.deleted_CDS_ID > ${i}.bin1.allchr.refExons.deleted
 echo ${i}.bin2.deleted_CDS_ID > ${i}.bin2.allchr.refExons.deleted
 echo ${i}.bin1.retained_CDS_ID > ${i}.bin1.allchr.refExons.retained
 echo ${i}.bin2.retained_CDS_ID > ${i}.bin2.allchr.refExons.retained
 for j in 01 02 03 04 05 06 07 08 09 10; do
 	cat Sb313_${i}_Chr${j}.bin1.noalignment.IDs.txt >> ${i}.bin1.allchr.refExons.noalign
 	cat Sb313_${i}_Chr${j}.bin2.noalignment.IDs.txt >> ${i}.bin2.allchr.refExons.noalign
 	cat Sb313_${i}_Chr${j}.bin1.refExons.deleted.IDs.txt >> ${i}.bin1.allchr.refExons.deleted
 	cat Sb313_${i}_Chr${j}.bin2.refExons.deleted.IDs.txt >> ${i}.bin2.allchr.refExons.deleted
 	cat Sb313_${i}_Chr${j}.bin1.refExons.retained.IDs.txt  >> ${i}.bin1.allchr.refExons.retained
 	cat Sb313_${i}_Chr${j}.bin2.refExons.retained.IDs.txt  >> ${i}.bin2.allchr.refExons.retained
 done
done

# QC check
for i in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	for b in bin1 bin2 ; do
		a=$(wc -l ${i}.${b}.allchr.refExons.noalign)
		d=$(wc -l ${i}.${b}.allchr.refExons.deleted)
		r=$(wc -l ${i}.${b}.allchr.refExons.retained)
		echo $i $b $a $d $r >> QC.totalcalls.txt
	done
done

sed -i 's/ /\t/g' QC.totalcalls.txt

#numbers come from wc -l ref_Sb313.cds.bed + 3 for the headers
#69327 + 3

awk -v OFS="\t" '{if ($3+$4+$5 == 69330) print $0, "TRUE" ; else print $0, "FALSE"}' QC.totalcalls.txt > QC.totalcalls.tsv
