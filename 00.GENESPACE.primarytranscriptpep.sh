#!/bin/bash

gff=$1
genome=$2
outname=$3

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/

ml singularity

#use this one if gff3 DOES specify canonical transcript in descripter field (field 9)
awk -F "\t" '$3=="mRNA" {print $9}' ${gff} | grep "canonical_transcript=1" | cut -f 1-2 -d ";" |sed 's/;Parent=/\t/g' |sed 's/ID=//g' |awk '{print $1"\t"$2}' > ${outname}-primary-transcript-ids.txt

#uses mikado utils grep to pull out primary transcripts without negating the features/nesting that makes a gff a gff
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep ${outname}-primary-transcript-ids.txt ${gff} ${outname}-primary-transcript.gff

ml samtools
ml cufflinks
samtools faidx ${genome}

gffread ${gff} -g ${genome} -x ${outname}.cds.fasta -y ${outname}.pep.fasta

ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' ${outname}.pep.fasta > ${outname}.pep.faa
