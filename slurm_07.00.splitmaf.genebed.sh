#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=1:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml bedops
bedops gff2starch < /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3
bedops gff2bed < /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 > /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.gene.bed

