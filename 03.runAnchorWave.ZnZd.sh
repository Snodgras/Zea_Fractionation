#!/bin/bash

FILE1=$1
FileName=${FILE1#*/}
GenomeName=$2
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \
	-as Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
	-a sam_files/${FileName}_Sb313.cds.sam \
	-ar sam_files/Sb313.ref.sam \
	-s ${FILE1}.fasta \
	-n AnchorWave_output/Sb313_${GenomeName}_anchorwave.anchors \
	-R 4 \
	-Q 1 \
	-o AnchorWave_output/Sb313_${GenomeName}_anchorwave.4to1.maf \
	-t 5 > AnchorWave_logs/Sb313_${GenomeName}_anchorwave.4to1.log 2>&1
