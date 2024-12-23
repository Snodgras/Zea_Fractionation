#!/bin/bash
ml singularity
genome=$1
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName swap_${genome}/ Sb313_${genome}_anchorwave.swap.maf

