#!/bin/bash
module load python/2.7.18-2ut3ogj
genome=$1
cat Sb313_${genome}_anchorwave.maf | python2 /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/anchorwave-maf-swap.py > Sb313_${genome}_anchorwave.swap.maf

