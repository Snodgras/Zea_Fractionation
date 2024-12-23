#!/bin/bash
outname=$1
assembly=$2
bed=$3

#module load bedtools2
bedtools getfasta -fo split_genome_assemblies/${outname} -fi ${assembly} -bed ${bed}

