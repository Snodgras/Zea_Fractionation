#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --time=6:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml singularity
singularity exec --bind $PWD /ptmp/arnstrm/phg_latest.sif /tassel-5-standalone/run_pipeline.pl -VCFMetricsPlugin \
        	-vcfDir /work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/GVCF/ \
        	-outFile /work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/QC/GVCFmetrics.tsv \
        	-endPlugin
