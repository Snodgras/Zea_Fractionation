#!/bin/bash

REFfasta=$1 #Sorghum genome fasta
REFname=$2 #query name without extensions

##Convert maf to gvcf (this takes about 15 or so minutes):

#MAFToGVCFPlugin <options>
#-referenceFasta <Reference Fasta> : Input Reference Fasta (required)
#-mafFile <Maf File> : Input MAF file.  Please note that this needs to be a MAF file with 2 samples.  The first will be assumed to be the Reference and the second will be the assembly. (required)
#-sampleName <Sample Name> : Sample Name to write to the GVCF file as the genome header or ID (required)
#-gvcfOutput <Gvcf Output> : Output GVCF file name (required)
#-fillGaps <true | false> : When true, if the maf file does not fully cover the reference genome any gaps in coverage will be #filled in with reference blocks. This is necessary if the resulting GVCFs are to be combined. (Default: false)

module load singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/elena/.conda/phg_latest.sif /tassel-5-standalone/run_pipeline.pl \
	-Xmx300g \
	-MAFToGVCFPlugin \
	-referenceFasta ${REFfasta} \
	-mafFile Sb313_${REFname}_anchorwave.maf \
	-sampleName Sb313_${REFname} \
	-gvcfOutput Sb313_${REFname}_anchorwave.gvcf \
	-fillGaps true > ${REFname}.out_gvcf.log 2>&1
	

##The command "-Xmx300g" demands that Java has enough RAM for the job to be run; in this case, 300g

