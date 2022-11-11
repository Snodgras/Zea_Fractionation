#!/bin/bash

filenames=$1 #list of file names like *.fasta.masked

for sample in $filenames
        do
                echo $sample
                #describer=$(echo ${sample} | sed 's/.fasta.masked//')
                #echo $describer

trf ${sample} 2 7 7 80 10 50 2000 -l 1 -d -h

done
