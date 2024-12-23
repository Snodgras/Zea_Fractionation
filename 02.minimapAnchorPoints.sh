#!/bin/bash
FILE1=$1 #Genomic fasta file ID for maize genome (genomeID)) (example: Zm-B73-REFERENCE-NAM-5.0); be sure to make certain the fasta file extension in the script matches user file extension i.e. '.fasta' vs '.fa', etc)
FileName=${FILE1#*/}

ml minimap2

#Do I need to redo the minimap to sorghum to itself everytime? Could cut down with a if file exists exemption
if [[ ! -f "sam_files/Sb313.ref.sam" ]] ; then
	minimap2 \
		-x splice \
		-t 11 \
		-k 12 \
		-a \
		-p 0.4 \
		-N 20 \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > sam_files/Sb313.ref.sam
fi

minimap2 \
	-x splice \
	-t 11 \
	-k 12 \
	-a \
	-p 0.4 \
	-N 20 \
	${FILE1} \
	Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > sam_files/${FileName%.fasta}_Sb313.cds.sam
