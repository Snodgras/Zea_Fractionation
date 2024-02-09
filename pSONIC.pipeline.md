#pSONIC pipeline

##0: reformat GFF
This step is to ensure we are only using the primary or canonical transcript for each genome's annotation. 
We also filter out any scaffolds that contain < 5 genes. 
This speeds up computational time in later steps. 
Because the Sorghum annotation was done by a different group, the gff file is formatted differently.
I am using the file `Sbicolor_313_v3.1.cds+primaryTranscriptOnly.fa.gz` from `https://data.jgi.doe.gov/refine-download/phytozome?q=313&expanded=Phytozome-313`

_Base Script_
`05.0.pSONIC.reformatGFF.sh`
```
gff=$1
genome=$2

ml singularity

#makes the id file for each gene's primary transcript 
#use this one if gff3 doesn't specify canonical transcript in descripter field (field 9)
#awk -F "\t" '$3=="mRNA" {print $9}' ${gff} | cut -f 1-2 -d ";" | sed 's/;Parent=/\t/g' |sed 's/ID=//g' |awk '{print $2"\t"$1}' |grep "_T001$" | awk '{print $2"\t"$1}' > ${genome}-primary-transcript-ids.txt

#use this one if gff3 DOES specify canonical transcript in descripter field (field 9)
awk -F "\t" '$3=="mRNA" {print $9}' ${gff} | grep "canonical_transcript=1" | cut -f 1-2 -d ";" |sed 's/;Parent=/\t/g' |sed 's/ID=//g' |awk '{print $1"\t"$2}' > ${genome}-primary-transcript-ids.txt

#uses mikado utils grep to pull out primary transcripts without negating the features/nesting that makes a gff a gff
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep ${genome}-primary-transcript-ids.txt ${gff} ${genome}-primary-transcript.gff

#filters out scaffolds  with less than 5 genes
## get scaffold names that have less than 5 genes
grep scaf ${genome}-primary-transcript.gff | awk -v OFS='\t' '$3=="gene" {print $1}' - | sort | uniq -c | sed 's/ //g' - | tr '[0-9]s' '[0-9]\t'| sed 's/caf/scaf/g' - | awk -v OFS='\t' '$1 < 5 {print $2}' - > ${genome}.scaf.ids
## use grep to exclude scaffolds that don't pass
grep -vwF -f ${genome}.scaf.ids ${genome}-primary-transcript.gff > ${genome}-primary-5andup.gff3

```

_Slurm Script_
_For 2 genomes, ~2.5 minutes_
_For 5 genomes, ~2.66 minutes_
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#unzip annotation files
#gunzip annotations_final/Zv-TIL01-REFERENCE-PanAnd-1.0_final.gff3.gz
gunzip annotations_final/Zx-TIL18-REFERENCE-PanAnd-1.0_final.gff3.gz
gunzip annotations_final/Zd-Gigi-REFERENCE-PanAnd-1.0_final.gff3.gz
gunzip annotations_final/Td-FL_9056069_6-DRAFT-PanAnd-1.0_final.gff3.gz
#gunzip NAM-annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz

#for PanAnd Directory
for i in Zv-TIL01 Zx-TIL18 Zd-Gigi Td-FL_9056069_6 ; do
	gunzip annotations_final/${i}*.gff3.gz
	bash scripts/05.0.pSONIC.reformatGFF.sh annotations_final/${i}-*.gff3 ${i}
done

#for NAM Directory
for i in Zm-B73 ; do
	gunzip NAM-annotations/${i}*.gff3.gz
	bash scripts/05.0.pSONIC.reformatGFF.sh NAM-annotations/${i}-*.gff3 ${i}
done

#for sorghum directory
#Sorghum 313 already has a primary transcript file associated with it, so we don't need to make it

gunzip Sbicolor_313_v3.1.cds_primaryTranscriptOnly.fa.gz
grep "^>" Sbicolor_313_v3.1.cds_primaryTranscriptOnly.fa | cut -f 4-5 -d " " | sed 's/locus=//g' - | sed 's/ID=//g' - | awk -v OFS='\t' '{print $1,$2.v3.1}' - > Sb-313-v3-primary-transcript-ids.txt

ml singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep Sb-313-v3-primary-transcript-ids.txt Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 Sb-313-v3-primary-transcript.gff

#filters out scaffolds  with less than 5 genes
## get scaffold names that have less than 5 genes
grep super Sb-313-v3-primary-transcript.gff | awk -v OFS='\t' '$3=="gene" {print $1}' - | sort | uniq -c | sed 's/ //g' - | tr '[0-9]s' '[0-9]\t'| sed 's/uper/super/g' - | awk -v OFS='\t' '$1 < 5 {print $2}' - > Sb-313-v3.scaf.ids

## use grep to exclude scaffolds that don't pass
grep -vwF -f Sb-313-v3.scaf.ids Sb-313-v3-primary-transcript.gff > Sb-313-v3-primary-5andup.gff3
```

*Quality Control Check*
Make sure that the number of genes pulled out of the original gff matches your expectations
```
grep -v "^#" {genome}-primary-5andup.gff3 |cut -f 3 |sort |uniq -c
```

I moved all of the files into `primarytranscriptfiles` and then moved all the `5andup.gff3` files back to the working directory to keep things neat. 

##1: Convert GFF to proteome fasta

_Base Script_
`05.1.pSONIC.orthofinderInputs.sh`
```
genome=$1
gff=$2
outname=$3

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/

#make a fasta file index for faster jobs
module load samtools
samtools faidx ${genome}

#get cds and protein fasta from gff and reference genome
module load cufflinks
gffread ${gff} -g ${genome} -x ${outname}.cds.fasta -y ${outname}.pep.fasta

```

_Slurm Script_
_~1.5 min_
`slurm_05.1.pSONIC.orthofinderInputs.sh`
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#for PanAnd Directory
for i in Zv-TIL01 Zx-TIL18 Zd-Gigi Td-FL_9056069_6 ; do
	gunzip assemblies_final/${i}*.fasta.gz
done

#for NAM Directory 
for i in Zm-B73 ; do
	gunzip NAM-annotations/${i}-*.gff3 ${i}
done

mkdir cds-fastas
mkdir pep-fastas

#for PanAnd Directory
for i in Zv-TIL01 Zx-TIL18 Zd-Gigi Td-FL_9056069_6 ; do
	bash scripts/05.1.pSONIC.orthofinderInputs.sh assemblies_final/${i}*.fasta ${i}-primary-5andup.gff3 ${i} ; 
done

#for NAM Directory
for i in Zm-B73 Zm-CML333 Zm-NC358 Zm-Oh43 ; do
	bash scripts/05.1.pSONIC.orthofinderInputs.sh NAM-assemblies/${i}*.fasta ${i}-primary-5andup.gff3 ${i} ; 
done

#for sorghum directory
bash scripts/05.1.pSONIC.orthofinderInputs.sh Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa Sb-313-v3-primary-5andup.gff3 Sb-313-v3

mv *cds.fasta cds-fastas/.
mv *pep.fasta pep-fastas/.

ml bioawk
cd pep-fastas
for i in *.pep.fasta ; 
	do bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' $i > ${i%.fasta}.faa ;  
done
```

*Quality Control Check*
Check the line count of the .faa files, should be 2x the number of genes in the gff
_I will need to change the bioawk step in the above script so that it'll run in a loop...won't right now_
_Probably due to the remove ${i} command at the end (removed that line so we'll see if it runs correctly now_

##2. Run Orthofinder

_Base Script_
`05.2.pSONIC.runOrthofinder.sh`
```
pepdir=$1

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
ml miniconda3

source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate orthofinder

orthofinder -og -t $SLURM_JOB_CPUS_PER_NODE -f ${pepdir} -S diamond
```

_Slurm Script_
_~6.5 minutes on 2 genomes_
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --exclusive   # take over the nodes processors and memory - all for yourself 
#SBATCH --partition=amd
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/05.2.pSONIC.runOrthofinder.sh pep-fastas/
```

*Quality Control Check*
Look at the slurm output file. It'll list how many orthogroups it made in the run plus other statistics.
This is a great chance to see if there was weirdness in running orthofinder.

##3. Prepare directory and files for MCScanX
###Part 1
_Base Script_
`05.3.1.pSONIC.makeCombinedGFF.sh`
```
rundate=$1
prefix=$2

mkdir ${prefix}_pSONIC_results
cp pep-fastas/OrthoFinder/${rundate}/WorkingDirectory/SpeciesIDs.txt ${prefix}_pSONIC_results/.
cp pep-fastas/OrthoFinder/${rundate}/WorkingDirectory/SequenceIDs.txt ${prefix}_pSONIC_results/.
cp pep-fastas/OrthoFinder/${rundate}/Orthogroups/Orthogroups.tsv ${prefix}_pSONIC_results/.
for i in pep-fastas/OrthoFinder/${rundate}/WorkingDirectory/Blast*  ; do cp $i ${prefix}_pSONIC_results/. ; done

cd ${prefix}_pSONIC_results
gunzip Blast*.gz
cat Blast*.txt > ${prefix}.blast
cd ..

#gets them into the basic order
for i in *5andup.gff3 ; do
	awk -v OFS='\t' '$3=="mRNA" {split($9,a,/;/); print $1, a[1], $4, $5}' ${i} > ${i}.temp
	sed -i 's/ID=//g' ${i}.temp
	grep -w scaf_[0-9] ${i}.temp
	grep -w scaf_10 ${i}.temp
	echo "Grepped for scaffolds numbered 1-10. If result above, need to rename."
done

```
*fasta has the "T001" suffix; must use mRNA*

_Slurm Script_
_~4 seconds_
`slurm_05.3.1.pSONIC.makeCombinedGFF.sh`
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --exclusive   # take over the nodes processors and memory - all for yourself 
#SBATCH --partition=amd
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/05.3.1.pSONIC.makeCombinedGFF.sh Results_Aug19 TIL01-TIL18

sed -i 's/chr/VA/g' Zv-TIL01-primary-5andup.gff3.temp
sed -i 's/scaf_/VA/g' Zv-TIL01-primary-5andup.gff3.temp
sed -i 's/chr/XA/g' Zx-TIL18-primary-5andup.gff3.temp
sed -i 's/scaf_/XA/g' Zx-TIL18-primary-5andup.gff3.temp

#Make sure this is going to ${prefix}.combined.gff in results directory
cat Zv-TIL01-primary-5andup.gff3.temp Zx-TIL18-primary-5andup.gff3.temp >> TIL01-TIL18_pSONIC_results/TIL01-TIL18.combined.gff

rm *.temp
```

###Part  2
_Base Script_
`05.3.2.pSONIC.translateGFF.sh`
```
prefix=$1
path2pSONIC=$2

ml miniconda3
source activate pSONICenv

cd ${prefix}_pSONIC_results
python ${path2pSONIC} ${prefix} translate_gff -gff ${prefix}.combined.gff
```

_Slurm Script_
_~3sec_
`slurm_05.3.2.pSONIC.translateGFF.sh`
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --exclusive   # take over the nodes processors and memory - all for yourself 
c#SBATCH --partition=amd
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/05.3.2.pSONIC.translateGFF.sh TIL01-TIL18 /work/LAS/mhufford-lab/snodgras/Fractionation/pSONIC/pSONIC.py
```
*Should make a file called <prefix>.gff*

##4. Run MCScanX

_Base Script_
`05.4.pSONIC.mcscanx.sh`
```
prefix=$1
path2MCScanX=$2

ml miniconda3
source activate pSONICenv

cd ${prefix}_pSONIC_results
${path2MCScanX} -b 2 ${prefix}
```

_Slurm Script_
_~34 sec_
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --exclusive   # take over the nodes processors and memory - all for yourself 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/05.4.pSONIC.mcscanx.sh TIL01-TIL18 /work/LAS/mhufford-lab/snodgras/Fractionation/MCScanX/MCScanX
```

##5. Run pSONIC
_Base Script_
`05.5.pSONIC.psonic.sh`
```
prefix=$1
path2pSONIC=$2

ml miniconda3
source activate pSONICenv

cd ${prefix}_pSONIC_results

python ${path2pSONIC} ${prefix}
```

_Slurm Script_
_~43 sec_
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --exclusive   # take over the nodes processors and memory - all for yourself 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/05.5.pSONIC.psonic.sh TIL01-TIL18 /work/LAS/mhufford-lab/snodgras/Fractionation/pSONIC/pSONIC.py
```