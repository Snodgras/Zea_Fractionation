#!/bin/bash

FILE1=$1
FileName=${FILE1#NAM-assemblies/}
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \
	-as Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
	-a ${FileName}_Sbicolor_313_v3.0_cds.sam \
	-ar Sbicolor_313_v3.0_ref.sam \
	-s ${FILE1}.fasta \
	-n Sbicolor_313_v3.0_${FileName}_anchorwave.anchors \
	-R 2 \
	-Q 1 \
	-o Sbicolor_313_v3.0_${FileName}_anchorwave.maf \
	-t 5 > Sbicolor_313_v3.0_${FileName}_anchorwave.log 2>&1 
