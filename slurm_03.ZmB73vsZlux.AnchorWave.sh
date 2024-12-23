#!/bin/bash
#Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=96:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="AnchorWaveAlignment"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive

ml minimap2
minimap2 \
	-x splice \
	-t 11 \
	-k 12 \
	-a \
	-p 0.4 \
	-N 20 \
	NAM-assemblies/Zm-B73.fasta \
	NAM-annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 > sam_files/ZmB73.ref.sam


minimap2 \
	-x splice \
	-t 11 \
	-k 12 \
	-a \
	-p 0.4 \
	-N 20 \
	assemblies_final/Zl-RIL003-REFERENCE-PanAnd-1.0.fasta  \
	NAM-assemblies/Zm-B73.fasta  > sam_files/ZlRIL003_ZmB73.cds.sam

ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i NAM-annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
	-as NAM-assemblies/Zm-B73.fasta \
	-r NAM-assemblies/Zm-B73.fasta \
	-a sam_files/ZlRIL003_ZmB73.cds.sam \
	-ar sam_files/ZmB73.ref.sam \
	-s assemblies_final/Zl-RIL003-REFERENCE-PanAnd-1.0.fasta \
	-n AnchorWave_output/ZmB73_ZlRIL003_anchorwave.anchors \
	-R 1 \
	-Q 1 \
	-o AnchorWave_output/ZmB73_ZlRIL003_anchorwave.maf \
	-t 5 > AnchorWave_logs/ZmB73_ZlRIL003_anchorwave.log 2>&1
