#!/bin/bash
FILE1=$1
FileName=${FILE1#NAM-assemblies/}
ml minimap2

minimap2 \
	-x splice \
	-t 11 \
	-k 12 \
	-a \
	-p 0.4 \
	-N 20 \
	${FILE1}.fasta \
	Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > ${FileName}_Sbicolor_313_v3.0_cds.sam

if [[ ! -f "Sbicolor_313_v3.0_ref.sam" ]] ; then

minimap2 -x splice -t 11 -k 12 -a -p 0.4 -N 20 Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > Sbicolor_313_v3.0_ref.sam

fi
