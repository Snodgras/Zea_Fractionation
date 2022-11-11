#pSONIC pipeline

##0: reformat GFF

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
	bash scripts/05.0.pSONIC.reformatGFF.sh annotations_final/${i}-*.gff3 ${i}
done

#for NAM Directory
for i in Zm-B73 ; do
	bash scripts/05.0.pSONIC.reformatGFF.sh NAM-annotations/${i}-*.gff3 ${i}
done

#add in unique modifiers to gene names so they're unique
#This step can be eliminated if the naming convention is unique for each genome
#sed -i 's/ID=/ID=Zv-TIL01_/g' Zv-TIL01-primary-5andup.gff
#sed -i 's/Parent=/Parent=Zv-TIL01_/g' Zv-TIL01-primary-5andup.gff
#sed -i 's/Name=/Name=Zv-TIL01_/g' Zv-TIL01-primary-5andup.gff

#sed -i 's/ID=/ID=Zx_TIL18_/g' Zx-TIL18-primary-5andup.gff
#sed -i 's/Parent=/Parent=Zx_TIL18_/g' Zx-TIL18-primary-5andup.gff
#sed -i 's/Name=/Name=Zx_TIL18_/g' Zx-TIL18-primary-5andup.gff
```

*Quality Control Check*
Make sure that the number of genes pulled out of the original gff matches your expectations
```
grep -v "^#" {genome}-primary-5andup.gff3 |cut -f 3 |sort |uniq -c
```

##1: Convert GFF to proteome fasta

_Base Script_
`05.1.pSONIC.orthofinderInputs.sh`
```
genome=$1
gff=$2
outname=$3

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

#for NAM Directory (WHICH NAM ASSEMBLIES TO USE MASKED OR UNMASKED?)
for i in Zm-B73 ; do
	gunzip NAM-annotations/${i}-*.gff3 ${i}
done

gunzip assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0.fasta.gz
gunzip assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta.gz

mkdir cds-fastas
mkdir pep-fastas

bash scripts/05.1.pSONIC.orthofinderInputs.sh assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0.fasta Zv-TIL01-primary-5andup.gff3 Zv-TIL01
bash scripts/05.1.pSONIC.orthofinderInputs.sh assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta Zx-TIL18-primary-5andup.gff3 Zx-TIL18

mv *cds.fasta cds-fastas/.
mv *pep.fasta pep-fastas/.

ml bioawk
cd pep-fastas
for i in *.pep.fasta ; 
	do bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' $i > ${i%.fasta}.faa ; 
	rm $i ;  
done
```

*Quality Control Check*
Check the line count of the .faa files, should be 2x the number of genes in the gff

##2. Run Orthofinder

_Base Script_
`05.2.pSONIC.runOrthofinder.sh`
```
pepdir=$1

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

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
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
#SBATCH --partition=amd
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