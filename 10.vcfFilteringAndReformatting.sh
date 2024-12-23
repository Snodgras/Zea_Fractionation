#!/bin/bash
vcf=$1 #like: chr10.bin2.headerless.indelonly.vcf 
output=${vcf%.headerless.indelonly.vcf}
chr=$2

#filter out insertions
head -n 1 ${vcf} > ${output}.delonly.vcf
cat ${vcf}| \
awk -v OFS='\t' '{split($5,a,/,/); if(length(a[1]) < length($4) || (length(a[2]) < length($4) && a[2] != "")) print $0}' - >> ${output}.delonly.vcf

#makes it a bed with the stop being the end of the REF allele
awk -v OFS='\t' '{split($5,a,/,/); 
		print $1,$2,$2+length($4),$4,$5,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45}' ${output}.delonly.vcf > ${output}.delonly.bed

ml bedtools2

#create exonic 
bedtools intersect -wa -wb -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${chr}.bed -b ${output}.delonly.bed > ${output}.del.exonic.bed

#create exact exonic
bedtools intersect -wa -wb -F 1 -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${chr}.bed -b ${output}.delonly.bed > ${output}.del.exact.exonic.bed

