# FRACTIONATION PIPELINE
# Maggie's Original Pipeline
**Key steps:**

1) Find tandem repeats in each masked genome using Tandem Repeat Finder
2) Align Sorghum exons associated with B73 subgenomes to each masked genome using dc-megablast
3) Remove dc-megablast alignments that overlap tandem repeat genomic coordinates as determined by Tandem Repeat Finder
4) Use the tandem-repeat filtered dc-megablast alignments from step 3) as inputs to DagChainer and run DagChainer
5) Organize and filter DagChainer outputs and generate tables


**Input files:**

a) Sbicolor_313_v3.1_exons_primary_notandems_cshl_clusters4.fa <br/>
b) Sb_exons_coords_CSHL_subgenomes_sections_sorted.txt<br/>
c) Repeatmasked NAM and B73 genomes<br/>


These steps were all done recursively on each masked NAM genome and B73.


**Software used for these analyses:**

a) Tandem Repeat Finder (TRF) 4.09 https://tandem.bu.edu/trf/trf.html<br/>
b) Blast+ 2.7.1<br/>
c) DagChainer http://dagchainer.sourceforge.net/<br/>
d) bedtools2 2.27 https://github.com/arq5x/bedtools2<br/>
e) R 3.6.0 (for post-analysis processing)<br/>
	- dplyr
	- purrr
	- tibble
	- tidyr
	- readr


## 1) Find Tandem Repeats


### Step 1: Run Tandem Repeat Finder (TRF)

```bash
sh trf.sh
```


**trf.sh** code

```bash
#!/bin/bash

for sample in *.fasta.masked
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked//')
                echo $describer #gets the name of the genome you're working on

trf ${sample} 2 7 7 80 10 50 2000 -l 1 -d -h

done
```
Arguments are:
match, mismatch, and delta weights: 2, 7, 7 (recommended values)
PM and PI: 80, 10 (recommended for best performance)
Minscore: 50 (example value)
Maxperiod: 2000 (max size the program can find)
-l 1: longest tandem array in million bp (default is 2, but may be too permissible for memory)
-d: produces a data file (text)
-h: suppresses HTML output

### Step 2: Parse TRF output:


```bash
sh parse_trf_bedfile_allsizes.sh
```	

**parse_trf_bedfile_allsizes.sh** code

```bash
#!/bin/bash

for sample in *.fasta.masked.2.7.7.80.10.50.2000.dat
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked.2.7.7.80.10.50.2000.dat//')
                echo $describer

#Convert trf .dat file to .bed file using trf's TRFdat_to_bed.py script:

python TRFdat_to_bed.py --dat ${sample} --bed ${describer}.bed

#Parse above bed file:

cut -f 1,2,3 ${describer}.bed > ${describer}_trf_allsizes.bed

done
```




## 2) Align Sorghum exons to each masked genome using dc-megablast

### Step 1: Create blast database:

```bash
sh blastdb.sh
```

**blastdb.sh** code

```bash
#!/bin/bash

module load blast-plus

#looping blastdb:

for sample in *.fasta.masked
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked//')
                echo $describer

makeblastdb -in ${sample} -dbtype nucl -out ${describer}_masked

done
```


### Step 2: align reads using dc-megablast:


```bash
sh dc-megablast.sh
```

**dc-megablast.sh** code

```bash
#!/bin/bash

module load blast-plus

# have all your masked blast databases within the same directory

for database in *.nhr

        do

blastn -task dc-megablast -outfmt 6 -query Sbicolor_313_v3.1_exons_primary_notandems_cshl_clusters4.fa -db ${database%.*} -num_threads 4 -out "${database}_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt"

        done
```


### Step 3: parse dc-megablast output:

```bash
sh parse_dc-megablast.sh
```

**parse_dc-megablast.sh** code

```bash
#!/bin/bash

#This script formats the dc-megablast outfmt 6 output and merges it with the coordinates associated with the Sorghum exons based on Sorghum exon ID for preparation for DagChainer; DagChainer needs a very specific format, as specified in their documentation. Sorghum exons not on chromosomes are filtered out. The Sorghum exons in Sb_exons_coords_CSHL_subgenomes_sections_sorted.txt have the subgenome and chromosome information vs B73 associated with them; alignments to the NAMs to Sorghum that don't match the B73 chr ID were ignored. 

module load bedtools2

for sample in *_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt//')
                echo $describer

cat ${sample} | grep -v "Sobic.K" | grep -v "super" | awk -v OFS="\t" '{if($9>$10) print $1,$2,$10,$9,$11; else print $1,$2,$9,$10,$11}' | sort -k1,1 | join - Sb_exons_coords_CSHL_subgenomes_sections_sorted.txt | tr ' ' '\t' | awk -v OFS="\t" '{if($2==$11) print $6,$1"_"$9"_"$10"_"$11,$7,$8,$2,$2"_"$3"_"$4,$3,$4,$5}' > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer.txt

	done
```




## 3) Remove dc-megablast alignments that overlap Tandem Repeat Finder coordinates

```bash
sh intersect_trf_allsizes.sh
```


**intersect_trf_allsizes.sh** code

```bash
#!/bin/bash

#This script compares the coordinates of the tandem array regions determined by trf with the coordinates of the filtered, formatted dc-megablast outputs and reports only non-intersecting dc-megablast alignment coordinates using bedtools intersect. The formatted dc-megablast output is first reformatted to a bed file, then formatted back to its original state after the bedtools step. This step could be optimized by formatting the dc-megablast output to a bedfile from the outset. 

module load bedtools2

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer.txt//')
                echo $describer


cat ${sample} | awk -v OFS="\t" '{print $5,$7,$8,$1,$2,$3,$4,$5,$9}' | bedtools sort | bedtools intersect -a - -b ${describer}*_trf_allsizes.bed -v | awk -v OFS="\t" '{print $4,$5,$6,$7,$1,$1"_"$2"_"$3,$2,$3,$9}' > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt

done
```



## 4) Use the tandem-repeat filtered dc-megablast alignments from step 3) as inputs to DagChainer and run DagChainer


```bash
sh dagchainer-filtered-trf.sh
```

**dagchainer-filtered-trf.sh** code

```bash
#!/bin/bash

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt//')
                echo $describer


#First, do another filter for tandem repeats using Dagchainer's filter_repetitive_matches.pl script:

perl /DAGCHAINER/accessory_scripts/filter_repetitive_matches.pl 100000 < ${sample} > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered;

#Run dagchainer (the output will append the extension .aligncoords to each output file); parameters were optimized for maize, and -A 15 takes into consideration the fact that individual exons, not genes, are being processed:

perl /DAGCHAINER/run_DAG_chainer.pl -i ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered -s -I -D 1000000 -g 40000 -A 15

        done
```



## 5) Organize and filter DagChainer outputs and generate table

This generates the exon count table of syntenic Sorghum aligned exons vs B73 and the NAMs within each syntenic block designated by DagChainer. Files were not sorted prior to processing since the syntenic block order from DagChainer was used. Exon counts were determined using bedtools groupby by grouping exons within a syntenic region that had the same Sorghum gene ID in that region and reporting the count of those exons, as well as the chromosome ID and the minimum coordinate (i.e. what results in the start coordinate) of those exons within the syntenic region. The result of this parsing was made into a table using R.


### Step 1: Get exon counts for each NAM genome and B73:

```bash
sh post_dagchainer_coords_csv_filtered_subg_unique_trf.sh
```

**post_dagchainer_coords_csv_filtered_subg_unique_trf.sh** code

```bash
#!/bin/bash

#Each subgenome was processed individually. Data was converted to csv files. 

module load bedtools2

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered.aligncoords
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered.aligncoords//')
                echo $describer

cat ${sample} | grep -v "#" | grep "M1" | awk '{$1=$1} 1' FS=".1.v3.1." OFS="\t" | tr '_' '\t' | awk -v OFS="\t" '{print $2".1_"$4"_"$5,$10,$11}' | groupBy -i - -g 1,2 -c 1,3 -o count,min  | sort -k1,1 -k3,3r | awk '!x[$1]++' | awk -v x=$describer 'BEGIN { OFS="\t"; print "gene_chr","pos_"x,"exoncount_"x} { print $1"_"$2,$4,$3}' | tr '\t' ',' > ${describer}_ISUmasked_Sb_dagchainer_filtered_dc-m_coords_trf_M1.csv

cat ${sample} | grep -v "#" | grep "M2" | awk '{$1=$1} 1' FS=".1.v3.1." OFS="\t" | tr '_' '\t' | awk -v OFS="\t" '{print $2".1_"$4"_"$5,$10,$11}' | groupBy -i - -g 1,2 -c 1,3 -o count,min  | sort -k1,1 -k3,3r | awk '!x[$1]++' | awk -v x=$describer 'BEGIN { OFS="\t"; print "gene_chr","pos_"x,"exoncount_"x} { print $1"_"$2,$4,$3}' | tr '\t' ',' > ${describer}_ISUmasked_Sb_dagchainer_filtered_dc-m_coords_trf_M2.csv

done
```

### Step 2: Make table for each subgenome:

all M1 csv files were put into their own directory; same with all M2 csv files.

For each directory:

Module load r

```bash
./r_table.R
```

**r_table.R** code

```bash
#!/usr/bin/env Rscript

#make a table from all csv files in a directory

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(readr)

multmerge = function(mypath){filenames = list.files(path=mypath, full.names=TRUE)
datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

full_data = multmerge("/path/to/directory")

write.csv(full_data, file = "output.csv")
```


### Step 3: Join subgenome tables:

Subgenome tables were joined together using unix join.

#########################################

##Running GENESPACE to identify orthologous blocks between/within Tripsacinae genomes
"https://github.com/jtlovell/GENESPACE"

Using the information from Arun's github (`https://github.com/HuffordLab-Containers/genespace_docker`)
Start at step 3 to put container on hpc
Add in the line to use the old modules for the slurm script in step 4, but otherwise copy exactly
`module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/`
Once the script is running on hpc from step 4, open browser to follow step 5 to work interactively through Rstudio

But first I need to set up the directory hierarchy and check the gff file formats to make sure it'll parse correctly
I'm comparing to the files Arun used to run Genespace:`/ptmp/arnstrm/To_Sam/GeneSpace_run1`

*Sorghum Genome Files*
This is where the sorghum genome and annotations come from: 
`https://data.jgi.doe.gov/refine-download/phytozome?organism=Sbicolor&expanded=Phytozome-313`

*Genome Names can only be alphanumeric characters, no hyphens*
Set up directories
```
/{RUN_NAME}/rawGenomes/{GENOME_NAME}/{GENOME_NAME}/annotation/.gff3
cd /GENESPACE_prelim/rawGenomes
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/Td-FL_9056069_6-DRAFT-PanAnd-1.0_final.gff3 Td-FL/Td-FL/annotation/
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/Zd-Gigi-REFERENCE-PanAnd-1.0_final.gff3 Zd-Gigi/Zd-Gigi/annotation/
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/Zv-TIL01-REFERENCE-PanAnd-1.0_final.gff3 Zv-TIL01/Zv-TIL01/annotation/
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/Zx-TIL18-REFERENCE-PanAnd-1.0_final.gff3 Zx-TIL18/Zx-TIL18/annotation/

ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 Zm-B73/Zm-B73/annotation
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-CML333-REFERENCE-NAM-1.0_Zm00026ab.1.gff3 Zm-CML333/Zm-CML333/annotation
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-NC358-REFERENCE-NAM-1.0_Zm00037ab.1.gff3 Zm-NC358/Zm-NC358/annotation
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/Zm-Oh43-REFERENCE-NAM-1.0_Zm00039ab.1.gff3 Zm-Oh43/Zm-Oh43/annotation
#DOUBLE CHECK WHICH GFF FOR SORGHUM TO USE
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.cds_primaryTranscriptOnly.fa Sb-313/Sb-313/annotation/
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 Sb-313/Sb-313/annotation/

```

Create primary transcript annotation files
_For PanAnd Genomes and NAM genomes_
Use `00.GENESPACE.primarytranscriptpep.sh` since they have the id tag for canonical transcript
For the slurm script use the commands:
```
cd GENESPACE_prelim/rawGenomes/
for i in Td-FL Zd-Gigi Zv-TIL01 Zx-TIL18 ; do 
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh ${i}/${i}/annotation/*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ${i}/${i}/annotation/${i} ;
done

for i in Zm-B73 Zm-CML333 Zm-NC358 Zm-Oh43 ; do 
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh ${i}/${i}/annotation/*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ${i}/${i}/annotation/${i} ;
done
```

`00.GENESPACE.primarytranscriptpep.sh`
```
#!/bin/bash

gff=$1
genome=$2
outname=$3

module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/

ml singularity

#use this one if gff3 DOES specify canonical transcript in descripter field (field 9)
awk -F "\t" '$3=="mRNA" {print $9}' ${gff} | grep "canonical_transcript=1" | cut -f 1-2 -d ";" |sed 's/;Parent=/\t/g' |sed 's/ID=//g' |awk '{print $1"\t"$2}' > ${outname}-primary-transcript-ids.txt

#uses mikado utils grep to pull out primary transcripts without negating the features/nesting that makes a gff a gff
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep ${outname}-primary-transcript-ids.txt ${gff} ${outname}-primary-transcript.gff

ml samtools
ml cufflinks
samtools faidx ${genome}

gffread ${gff} -g ${genome} -x ${outname}.cds.fasta -y ${outname}.pep.fasta

ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' ${outname}.pep.fasta > ${outname}.pep.faa
```
_For Sorghum genome_
```
grep "^>" Sbicolor_313_v3.1.cds_primaryTranscriptOnly.fa | cut -f 4-5 -d " " | sed 's/locus=//g' - | sed 's/ID=//g' - | awk -v OFS='\t' '{print $2,$1}' - > Sb-313-v3-primary-transcript-ids.txt

ml singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep Sb-313-v3-primary-transcript-ids.txt Sbicolor_313_v3.1.gene.gff3 Sb-313-v3-primary-transcript.gff

ml samtools
ml cufflinks
samtools faidx /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa
gffread Sb-313-v3-primary-transcript.gff -g /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa -x Sb-313.cds.fasta -y Sb-313.pep.fasta
ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' Sb-313.pep.fasta > Sb-313.pep.faa
```
_Clean pep file IDs_
Need to look at sorghum
Arun has a script called `cleanPep.sh` in his `rawGenomes` directory
I need to remove the `_T001` from the ends of the PanAnd and NAM fasta IDs
I need to change sorghum IDs from `Sobic.001....1.v3.1` to `Sobic001G...` (no periods)
```
#!/bin/bash
file=$1
if [ ${file} != Sb-313.pep.faa ]
then
	cut -f 1 -d "_" ${file} > temp
	mv temp ${file}
else
	cut -f 1-2 -d "." ${file} > temp
	sed -i 's/\.//g' temp
	mv temp ${file}
fi
```
To run the above script
```
for i in Sb-313 Td-FL Zd-Gigi Zm-B73 Zm-CML333 Zm-NC358 Zm-Oh43 Zv-TIL01 Zx-TIL18 ; do 
	cd ${i}/${i}/annotation/
	bash ../../../cleanpep.sh ${i}.pep.faa
	cd ../../..
done
```
Move all of the files except the `.pep` and `.gff` to a new directory called "extra" or genespace won't run
```
for i in T* Z* ; do cd ${i}/${i}/annotation ; mkdir extra ; mv *.fasta extra/. ; mv *.gff3 extra/. ; mv *.txt extra/. ; ls ; cd ../../.. ; done
```
_Start the interactive session_ (`sbatch scripts/genespace_container_v2.slurm`)
```
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:nova21-33:43619 snodgras@nova.its.iastate.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: snodgras
   password: qSEWuDQXYFFUvjIXOf37

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f 4389145

```
Errors thrown:
- parameters: Genome IDs have to be alphanumeric (no hyphens)
- parameters: Got warning messages about !grepl() having mismatching lengths and coercing to a logical(1)
- parse annotations: the gene IDs for sorghum don't match between the gff and the pep because of the pep.cleaning step
	Fixed with going back to the original pep.fasta and redoing the bioawk step
	then `cut -f 1-2,4-5 -d "." Sb-313.pep.faa > temp ; mv temp Sb-313.pep.faa`
- parse annotations: At the end of each parsing, there was the warning for grep and cat that there was a write error: Broken pipe
(Should still be running even if I close my browser and cut the connection via ssh, will have to check that later if the orthofinder step is still running)
- Note: orthofinder took ~30 minutes to run, MCScanX took (~30 minutes)
- after MCScanX: In any(genomeIDs) : coercing argument of type 'character' to logical

Output from the default run of the pangenome()

```
> pangenome(gsParam = gpar, refGenome = "Sb313")
*NOTE* RefGenome Sb313 has >1x ploidy - this is fine, but will be slower.
	Depending on the size of your run, you may run into memory issues.
Building reference-anchored scaffold against Sb313
	n. ref positions = 30234
	Reading in hits against Sb313 ... found 244805
	Interpolating positions ... n. genes mapped: 1x = 248290, 2+x = 32026, 0x = 47344
	Forming ref.-anchored db ... found 211951 genes for 33510 placements
Completing the pan-genome annotation ...
	Adding non-anchor entries ... found 17985 genes and 5167 placements
	Checking missing direct ref. syn. OGs ... found 3952 genes and 1632 placements
	Adding indirect syn. OGs ... found 492 genes and 240 placements
	Adding syn. OGs without ref. anchor ... found 76182 genes and 55640 placements
	Adding missing genes by synOG identity ... found 17098 genes and 376 placements
Annotating and formatting pan-genome
	Adding non-anchor entries ... found 24740 genes and 7755 placements
	No orthologue file available. Will ignore
	Writing pangenome to results/Sb313_pangenomeDB.txt.gz
	Returning wide-format with only syntenic array reps
	Done!
        pgID     pgChr   pgOrd                 Sb313            TdFL          ZdGigi           ZmB73        ZmCML333
    1:     1     Chr01     1.0 Sobic.001G000100.v3.1 Td00001aa023931 Zd00001aa005636 Zm00001eb065400 Zm00026ab065460
    2:     2     Chr01     2.0 Sobic.001G000200.v3.1 Td00001aa023932 Zd00001aa005635 Zm00001eb065390 Zm00026ab065450
    3:     3     Chr01     3.0 Sobic.001G000300.v3.1                                                                
    4:     4     Chr01     4.0 Sobic.001G000400.v3.1 Td00001aa023933 Zd00001aa005634 Zm00001eb065380 Zm00026ab065440
    5:     5     Chr01     4.5                                                       Zm00001eb065370                
   ---                                                                                                              
89145: 89145  super_88 33880.0    Sobic.K043600.v3.1                                                                
89146: 89146  super_88 33881.0    Sobic.K043700.v3.1                                                                
89147: 89147 super_896 33927.0    Sobic.K043800.v3.1                                                                
89148: 89148  super_99 33882.0    Sobic.K044200.v3.1                                                                
89149: 89149  super_99 33883.0    Sobic.K044300.v3.1                                                                
               ZmNC358          ZmOh43         ZvTIL01         ZxTIL18
    1: Zm00037ab065610 Zm00039ab065730 Zv00001aa005624 Zx00002aa005252
    2: Zm00037ab065600 Zm00039ab065720 Zv00001aa005623 Zx00002aa005251
    3:                                                                
    4: Zm00037ab065590 Zm00039ab065700 Zv00001aa005622 Zx00002aa005250
    5:                                                                
   ---                                                                
89145:                                                                
89146:                                                                
89147:                                                                
89148:                                                                
89149:  
```

###The 1.1.4 update in genespace

For each genome, I need: 
1. bed formatted coordinates of each gene (chr, start, end, name)
2. peptide sequences in fasta format, header exactly matches the "name" (4th) bed column

We still need to use the primary transcript, so run the `00.GENESPACE.primarytranscriptpep.sh`

Can use the partially automated annotation parsing (3.3 in the documentation)
Need the gff3 of protein coding gene features and the primary  transcript translated cds (peptides)
directory structure is
```
/genomeRepo
|__species1_genoX_v1.0_NCBI
    └─ peptidesGenoX.fa
    └─ genesGenoX.gff3
└─ species2_genoY_v1.0_NCBI
    └─ peptidesSpecies2.fa
    └─ genes.gff3
└─ species3_genoW_v1.0_phytozome
    └─ peptides.fa
    └─ genesSpecies2v1.gff3
└─ species4_genoZ_v1.0_otherRepo
    └─ peptides.fa
    └─ genes.gff3
  ...
```
Making file directory
```
mkdir GENESPACE_prelim.v1.1.4 
mkdir genomeRepo
for i in TdFL ZdGigi ZmB73 ZmCML333 ZmNC358 ZmOh43 ZvTIL01 ZxTIL18 Sb313 ; do mkdir genomeRepo/${i} ;done
#after running the 00.GENESPACE.primarytranscriptpep.sh in a separate directory 
for i in TdFL ZdGigi ZmB73 ZmCML333 ZmNC358 ZmOh43 ZvTIL01 ZxTIL18 Sb313 ; do
	cp GENESPACE_prelim/rawGenomes/${i}/${i}/annotation/*.pep.faa GENESPACE_prelim.v1.1.4/genomeRepo/${i}/.
	cp GENESPACE_prelim/rawGenomes/${i}/${i}/annotation/*primary-transcript.gff GENESPACE_prelim.v1.1.4/genomeRepo/${i}/.
done
```
Console output after doing the `run_genespace()`
```
GENESPACE run complete!  All results are stored in
/work/LAS/mhufford-lab/snodgras/Fractionation/GENESPACE_prelim.v1.1.4/ in the following subdirectories:
	syntenic block dotplots: /dotplots (...synHits.pdf)
	annotated blast files  : /syntenicHits
	annotated/combined bed : /results/combBed.txt
	syntenic block coords. : /results/blkCoords.txt
	syn. blk. by ref genome: /riparian/refPhasedBlkCoords.txt
	pan-genome annotations : /pangenes (...pangenes.txt.gz)
	riparian plots         : /riparian
	genespace param. list  : /results/gsParams.rda
```

### Prelim run 3
* Use everything as haploid assemblies like in the Genespace manuscript
* Use the psuedomolecule version of the Td genome
* Add in the Andropogon virginicus assembly as a good outgroup (potentially to use in place of Sorghum)

Making file directory
```
mkdir prelim.3
cd prelim.3
mkdir genomeRepo
for i in Td-FL Zd-Gigi Zm-B73 Zm-CML333 Zm-NC358 Zm-Oh43 Zv-TIL01 Zx-TIL18 Sb313 Av; do mkdir genomeRepo/${i} ;done
```
Run `00.GENESPACE.primarytranscriptpep.sh` on genomes from scratch using the following slurm script:
```
cd prelim.3/genomeRepo

#for pan-and genomes
for i in Td-FL Zd-Gigi Zv-TIL01 Zx-TIL18 Av; do
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3 ] 
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3.gz 
	fi
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ]
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta.gz
	fi	
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ${i}/${i};
done

#for NAM genomes
for i in Zm-B73 Zm-CML333 Zm-NC358 Zm-Oh43 ; do
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3 ] 
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3.gz
	fi
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ]
        then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta.gz
        fi
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ${i}/${i} ;
done

```
Then we need to make sure that in each genome directory, there's only 1 peptide fasta and 1 gff3 file
These are the files `*.pep.faa` and `*-primary-transcript.gff`

```
cd prelim.3/genomeRepo/
for i in * ; do
	rm ${i}/*.fasta
	rm ${i}/*ids.txt
done
```
The _sorghum genome_ denotes canonical transcript differently from the Pan-And and NAM genomes.
So use this to create the necessary files:
```
grep "^>" Sbicolor_313_v3.1.cds_primaryTranscriptOnly.fa | cut -f 4-5 -d " " | sed 's/locus=//g' - | sed 's/ID=//g' - | awk -v OFS='\t' '{print $2,$1}' - > Sb-313-v3-primary-transcript-ids.txt

ml singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep Sb-313-v3-primary-transcript-ids.txt Sbicolor_313_v3.1.gene.gff3 Sb-313-v3-primary-transcript.gff

ml samtools
ml cufflinks
samtools faidx /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa
gffread /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb-313-v3-primary-transcript.gff -g /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa -y Sb313/Sb313.pep.fasta
ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' Sb313/Sb313.pep.fasta > Sb313/Sb313.pep.faa

cp /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb-313-v3-primary-transcript.gff Sb313/Sb313-primary-transcript.gff

```
I need to remove the `_T001` from the ends of the PanAnd and NAM fasta IDs
I need to change sorghum IDs from `Sobic.001....1.v3.1` to `Sobic001G...` (no periods)
```
#!/bin/bash
file=$1
if [ ${file} != Sb313.pep.faa ]
then
	cut -f 1 -d "_" ${file} > temp
	mv temp ${file}
fi
```
To run the above script
```
for i in Av Td-FL Zd-Gigi Zm-B73 Zm-CML333 Zm-NC358 Zm-Oh43 Zv-TIL01 Zx-TIL18 ; do 
	cd ${i}/
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/cleanPep.sh ${i}.pep.faa
	cd ../
done
```

Note the genome names cannot have any non-alphanumeric characters aside form `_` and `.`

_Start the interactive session_ (`sbatch scripts/genespace_container_v2.slurm`)
```
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:nova21-33:43619 snodgras@nova.its.iastate.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: snodgras
   password: qSEWuDQXYFFUvjIXOf37

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f 4624463

```

I wrote the filtered dataframes as 
```
test.results.tsv #these are the coordinates that were found reciprocally exactly 1 time for a given block ID
#n=3401 unique sets of coordinates
test.lower.tsv #these are the genomes that did not have a reciprocal match or any match for a given block ID
#n=5456 unique sets of comparisons
test.higher.tsv #these are the genomes/coordinates that had several reciprocal matches for a given block ID
#n=368 sets of unique comparisons by block ID
```
_NEXT STEP_: Make bed file from the `test.results.tsv`? 


###Running the full set of genomes
* Include all NAM genomes
* Include all teosinte genomes

1. Make the directory
```
mkdir fullrun.1
cd fullrun.1
mkdir genomeRepo
for i in TdFL ZdGigi ZdMomo ZhRIMHU001 ZlRIL003 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 Sb313 Av ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ; do mkdir genomeRepo/${i} ;done
```
There should be 38 genome directories in the `/genomeRepo`

2. Make the primary peptide fastas
Run `00.GENESPACE.primarytranscriptpep.sh` on genomes from scratch using the following slurm script:
```
cd fullrun.1/genomeRepo

#for pan-and genomes
for i in Td-FL Zd-Gigi Zd-Momo Zh-RIMHU001 Zl-RIL003 Zn-PI615697 Zv-TIL01 Zv-TIL11 Zx-TIL18 Zx-TIL25 Av; do
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3 ] 
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3.gz 
	fi
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ]
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta.gz
	fi
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ${i//\-}/${i//\-};
done

#for NAM genomes
for i in Zm-B73 Zm-B97 Zm-CML103 Zm-CML228 Zm-CML247 Zm-CML277 Zm-CML322 Zm-CML333 Zm-CML52 Zm-CML69 Zm-HP301 Zm-IL14H Zm-Ki11 Zm-Ki3 Zm-Ky21 Zm-M162W Zm-M37W Zm-Mo18W Zm-MS71 Zm-NC350 Zm-NC358 Zm-Oh43 Zm-Oh7b Zm-P39 Zm-Tx303 Zm-Tzi8 ; do
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3 ] 
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3.gz
	fi
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ]
        then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta.gz
        fi
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/${i}*.fasta ${i//\-}/${i//\-} ;
done

```
Then we need to make sure that in each genome directory, there's only 1 peptide fasta and 1 gff3 file
These are the files `*.pep.faa` and `*-primary-transcript.gff`

```
cd fullrun.1/genomeRepo/
for i in * ; do
	rm ${i}/*.fasta
	rm ${i}/*ids.txt
done
```
The _sorghum genome_ denotes canonical transcript differently from the Pan-And and NAM genomes.
So use this to create the necessary files (if not already done so with the preliminary runs):
```
grep "^>" Sbicolor_313_v3.1.cds_primaryTranscriptOnly.fa | cut -f 4-5 -d " " | sed 's/locus=//g' - | sed 's/ID=//g' - | awk -v OFS='\t' '{print $2,$1}' - > Sb-313-v3-primary-transcript-ids.txt

ml singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/arnstrm/PanAnd/triffid/current-versions-genomes.triffid/mikado-v2.3.sif mikado util grep Sb-313-v3-primary-transcript-ids.txt Sbicolor_313_v3.1.gene.gff3 Sb-313-v3-primary-transcript.gff

ml samtools
ml cufflinks
samtools faidx /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa
gffread /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb-313-v3-primary-transcript.gff -g /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa -y Sb313/Sb313.pep.fasta
ml bioawk
bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' Sb313/Sb313.pep.fasta > Sb313/Sb313.pep.faa

cp /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb-313-v3-primary-transcript.gff Sb313/Sb313-primary-transcript.gff

```
I need to remove the `_T001` from the ends of the PanAnd and NAM fasta IDs
_DO NOT CHANGE THE SORGHUM GENE IDS OTHERWISE IT WONT WORK_
```
#!/bin/bash
file=$1
if [ ${file} != Sb313.pep.faa ]
then
	cut -f 1 -d "_" ${file} > temp
	mv temp ${file}
fi
```
To run the above script
```
for i in * ; do 
	cd ${i}/
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/cleanPep.sh ${i}.pep.faa
	cd ../
done
```

Note the genome names cannot have any non-alphanumeric characters aside form `_` and `.`

I had an issue with the parsing annotation step with the following genomes: 
ZhRIMHU001
ZlRIL003
ZmIL14H
ZmMS71
ZmOh7b

In looking at ZhRIMHU001 more closely, I noticed that there were duplicate transcript IDs in the pep file with different pep sequences
In doing the `00.GENESPACE.primarytranscriptpep.sh` step by step manually, the issue disappeared. 
I'm not sure why duplicates were present in the automated slurm script but not in doing it step by step. 

I'll do it step by step for the above genomes, using the following to double check for duplicates

```
wc -l ZhRIMHU001-primary-transcript-ids.txt 
uniq ZhRIMHU001-primary-transcript-ids.txt | wc -l
```

I'll replace the files in the `/genomeRepo` with the ones I create manually assuming they pass the check
Create an extra directory in their respective repos just in case we need the intermediate files

So the Zl genome failed because there was no gff file in `annotations_final` because we didn't generate RNA seq and thus have no annotations for that genome ::sad face::
Removing Zlux from the list

The others failed from NAM because the letters have different capitalizations!
Editing the slurm script to only redo those and changing the letter captialization 

3. Run GENESPACE

_Start the interactive session_ (`sbatch scripts/genespace_container_v2.slurm`)

```
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:nova21-32:46081 snodgras@nova.its.iastate.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: snodgras
   password: EMDJE6fHXd1Snbhxnpq9

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f 4648232

```

**NOTE** the genespace parameter object is returned or can be loaded into R via
        `load('/work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/results/gsParams.rda', verbose = TRUE)`. Then
        you can customize your riparian plots by calling `plot_riparian(gsParam = gsParam, ...)`. The source data and
        ggplot2 objects are also stored in the /riparian directory and can also be accessed by `load(...)`.
**NOTE** To query genespace results by position or gene, use `query_genespace(...)`. See specifications in
        ?query_genespace for details.
        
```
8. Constructing syntenic pan-gene sets ...
	Sb313     : n pos. = 170684, synOgs = 1246749, array mem. = 237475, NS orthos 1565850
	TdFL      : n pos. = 190475, synOgs = 1395255, array mem. = 264987, NS orthos 1758148
	ZdGigi    : n pos. = 187816, synOgs = 1403976, array mem. = 268162, NS orthos 1709021
	ZdMomo    : n pos. = 187426, synOgs = 1403235, array mem. = 268304, NS orthos 1707692
	ZhRIMHU001: n pos. = 187508, synOgs = 1400666, array mem. = 268374, NS orthos 1708448
	ZnPI615697: n pos. = 183290, synOgs = 1370676, array mem. = 262253, NS orthos 1670141
	ZvTIL01   : n pos. = 187704, synOgs = 1403520, array mem. = 268165, NS orthos 1708759
	ZvTIL11   : n pos. = 187553, synOgs = 1398421, array mem. = 267217, NS orthos 1707180
	ZxTIL18   : n pos. = 187460, synOgs = 1400616, array mem. = 268140, NS orthos 1708660
	ZxTIL25   : n pos. = 187552, synOgs = 1400390, array mem. = 267731, NS orthos 1707666
	ZmB73     : n pos. = 187897, synOgs = 1405288, array mem. = 269548, NS orthos 1703941
	ZmB97     : n pos. = 187843, synOgs = 1404997, array mem. = 270928, NS orthos 1700613
	ZmCML103  : n pos. = 187650, synOgs = 1401427, array mem. = 269660, NS orthos 1700665
	ZmCML228  : n pos. = 190239, synOgs = 1418360, array mem. = 272685, NS orthos 1720104
	ZmCML247  : n pos. = 187592, synOgs = 1398960, array mem. = 267431, NS orthos 1699367
	ZmCML277  : n pos. = 187767, synOgs = 1403564, array mem. = 270238, NS orthos 1700969
	ZmCML322  : n pos. = 187606, synOgs = 1403391, array mem. = 268375, NS orthos 1698569
	ZmCML333  : n pos. = 187723, synOgs = 1406851, array mem. = 270397, NS orthos 1696076
	ZmCML52   : n pos. = 188813, synOgs = 1411853, array mem. = 271140, NS orthos 1710974
	ZmCML69   : n pos. = 187782, synOgs = 1409120, array mem. = 270222, NS orthos 1699514
	ZmHP301   : n pos. = 187899, synOgs = 1402949, array mem. = 269295, NS orthos 1701443
	ZmIL14H   : n pos. = 187526, synOgs = 1401888, array mem. = 269297, NS orthos 1700165
	ZmKi11    : n pos. = 187556, synOgs = 1402894, array mem. = 269510, NS orthos 1697840
	ZmKi3     : n pos. = 187789, synOgs = 1401175, array mem. = 269084, NS orthos 1700850
	ZmKy21    : n pos. = 187772, synOgs = 1405764, array mem. = 269353, NS orthos 1698433
	ZmM162W   : n pos. = 187663, synOgs = 1398844, array mem. = 268590, NS orthos 1698587
	ZmM37W    : n pos. = 187603, synOgs = 1405140, array mem. = 270092, NS orthos 1699278
	ZmMo18W   : n pos. = 187836, synOgs = 1404564, array mem. = 269189, NS orthos 1699864
	ZmMS71    : n pos. = 188602, synOgs = 1410705, array mem. = 270442, NS orthos 1708919
	ZmNC350   : n pos. = 187693, synOgs = 1405864, array mem. = 269755, NS orthos 1697491
	ZmNC358   : n pos. = 187769, synOgs = 1404958, array mem. = 270482, NS orthos 1698787
	ZmOh43    : n pos. = 187672, synOgs = 1405095, array mem. = 271073, NS orthos 1699413
	ZmOh7b    : n pos. = 187608, synOgs = 1401352, array mem. = 268817, NS orthos 1698342
	ZmP39     : n pos. = 187485, synOgs = 1402411, array mem. = 269443, NS orthos 1697805
	ZmTx303   : n pos. = 188356, synOgs = 1406827, array mem. = 270383, NS orthos 1704587
	ZmTzi8    : n pos. = 188655, synOgs = 1407092, array mem. = 269684, NS orthos 1711445
	Av        : n pos. = 170836, synOgs = 1249239, array mem. = 238722, NS orthos 1568162

```

4. Split the genomes
Done also in R
Had to update the genespace image to 1.2.3 so it'd have the `plot_2genomes()`
But then needed minimap2 so updated the container again to `genespace_1.3.1.sif`
```
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:nova21-31:42095 snodgras@nova.its.iastate.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: snodgras
   password: pf11/wzIc8zkTG9F4oaa

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f 4654020
```

I need to create softlinks to the:
Genome Assemblies: /fafiles/[name].fasta
Gene Annotations: /gff3/[name].gff3
Repeat annotations: /repeats/[name].repeats.gff3

So they are named the same way for looping through 2genomes function

```
#make empty directories
cd fullrun.1
mkdir fafiles
mkdir gff3
mkdir repeats

for i in /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-assemblies/*.fasta ; do ln -s ${i} fafiles/. ; done
cd fafiles
for i in *.fasta ; do mv ${i} Zm${i:3} ; done
for i in /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/*.fasta ; do ln -s ${i} . ; done
rm Zl-RIL003-REFERENCE-PanAnd-1.0.fasta
for i in *REFERENCE*fasta ; do mv ${i} ${i:0:2}${i:3:8}.fasta ;done
mv AvKellogg1 Av.fasta
mv TdFL_90560 TdFL.fasta
for i in *REF ; do mv ${i} ${i:0:6}.fasta ; done
for i in *RE ; do mv ${i} ${i:0:7}.fasta ; done
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa Sb313.fasta

cd ../gff3
for i in /work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/* ; do s=$(echo ${i#/work/LAS/mhufford-lab/snodgras/Fractionation/NAM-annotations/}) ; t=$(echo ${s%-REFERENCE*}.gff3) ; ln -s ${i} Zm${t:3} ; done
mv ZmIl14H.gff3 ZmIL14H.gff3
mv ZmMs71.gff3 ZmMS71.gff3
mv ZmOh7B.gff3 ZmOh7b.gff3
for i in /work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/*.gff3 ; do s=$(echo ${i#/work/LAS/mhufford-lab/snodgras/Fractionation/annotations_final/} ); ln -s ${i} ${s%-REFERENCE*}.gff3 ; done
mv Av-Kellogg1287_8.gff3 Av.gff3
mv Td-FL_9056069_6.gff3 TdFL.gff3
for i in Zd* ; do mv ${i} Zd${i:3}; done
for i in Zv* Zx* ; do mv ${i} ${i:0:2}${i:3} ;done
for i in *-* ; do mv ${i} ${i:0:2}${i:3} ;done
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 Sb313.gff3 

for i in /work/LAS/mhufford-lab/snodgras/Fractionation/repeats_final/* ; do ln -s $i . ; done
mv Sbicolor_313_v3.1.repeatmasked_assembly_v3.0.gff3 Sb313.repeats.gff3
for i in *REFERENCE* ; do s=$(echo ${i%-REFERENCE*}.repeats.gff3) ;mv ${i} ${s:0:2}${s:3} ;done
mv AvKellogg1287_8.repeats.gff3 Av.repeats.gff3
mv TdFL_9056069_6.repeats.gff3 TdFL.repeats.gff3
mv ZmIl14H.repeats.gff3 ZmIL14H.repeats.gff3
mv ZmOh7B.repeats.gff3 ZmOh7b.repeats.gff3

```

Trouble shooting why it's not working with ZdGigi or Zea nica.
```
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/fullrun.1/genomeRepo/ZdGigi/ZdGigi-primary-transcript.gff ZdGigi.primarytranscript.gff3
#soft link all the assemblies and repeat annotations to also be the same root name

#try removing all the scaffolds and running with just chromosomes
for i in {1..10}; do echo chr$i; done > ids-to-keep.txt
seqtk subseq ZdGigi.fasta ids-to-keep.txt > filtered.ZdGigi.fasta
grep -Ev "^scaf|^alt-scaf" ZdGigi.gff3 > filtered.ZdGigi.gff3
#For Zn, change alt-scaf to alt-ctg and scaf to ctg
grep -Ev "^ctg|^alt-ctg" ZnPI615697.gff3
```

4. Figure out the coordinates for the duplicate regions and split by synteny

Use the `*.blockCoords.txt` from the `plot_2genomes()`
What it looks like:
```
(base) [snodgras@nova output]$ head Sb313_vs_ZmB73.blockCoords.txt 
blkID	genome1	genome2	chr1	chr2	start1	end1	min2	max2	nHits	orient	start2	end2
Chr03_chr3_3	Sb313	ZmB73	Chr03	chr3	25240	1813144	1212517	3612786	481	+	1212517	3612786
Chr10_chr9_12	Sb313	ZmB73	Chr10	chr9	13192	91466	30849476	31479988	43	+	30849476	31479988
Chr02_chr7_1	Sb313	ZmB73	Chr02	chr7	138706	1071539	180823	2340658	283	+	180823	2340658
Chr05_chr4_2	Sb313	ZmB73	Chr05	chr4	89730	285193	188383558	188701083	35	+	188383558	188701083
Chr04_chr5_6	Sb313	ZmB73	Chr04	chr5	117387	134487	68402642	68421148	32	+	68402642	68421148
Chr10_chr9_1	Sb313	ZmB73	Chr10	chr9	99209	115263	32132232	32173221	28	+	32132232	32173221
Chr04_chr5_29	Sb313	ZmB73	Chr04	chr5	331163	537213	69155462	69798951	62	+	69155462	69798951
Chr08_chr10_2	Sb313	ZmB73	Chr08	chr10	316333	2051296	581056	2523428	230	+	581056	2523428
Chr05_chr4_17	Sb313	ZmB73	Chr05	chr4	355469	743586	188949293	189299895	52	+	188949293	189299895
```

Going to test the script using just chr10 of Sb313:
```
####STEP 1.1: Split query genome (ZmB73) by chromosomes 
grep -E "Chr10" Sb313_vs_ZmB73.blockCoords.txt > test.blockCoords.txt
#works with manually specifying
if grep -Eq "chr1" test.blockCoords.txt 
then
	grep -E "chr1" test.blockCoords.txt > chr1.test.blockCoords.txt
fi
#get it to work in a loop
for i in {1..10} ; do
	if grep -Eq $(echo chr${i}) test.blockCoords.txt
	then
		grep -E $(echo chr${i}) test.blockCoords.txt > chr${i}.test.blockCoords.txt
	fi 
done

####Step 1.2: convert to bed file format:
for i in 5 6 9 ; do
awk -v OFS='\t' '{print $4,$6,$7,$3,$5,$12,$13}' chr${i}.test.blockCoords.txt > chr${i}.test.bed
done

####STEP1.3: Intersect chromosome files
module load bedtools2
bedtools intersect -wa -wb -a chr6.test.bed -b chr9.test.bed >> intersect.test.bed
bedtools intersect -wa -wb -a chr5.test.bed -b chr9.test.bed >> intersect.test.bed
bedtools intersect -wa -wb -a chr6.test.bed -b chr5.test.bed >> intersect.test.bed
```
Ok, it was weird with a lot of Sb313 chr10 only having 1 syntenic match (missing duplicate).

Testing on a different chromosome:
```
####STEP 1.1: Split query genome (ZmB73) by chromosomes 
grep -E "Chr02" Sb313_vs_ZmB73.blockCoords.txt > chr2test.blockCoords.txt
for i in {1..10} ; do
	if grep -Eq $(echo chr${i}) chr2test.blockCoords.txt
	then
		grep -E $(echo chr${i}) chr2test.blockCoords.txt > chr${i}.chr2test.blockCoords.txt
	fi 
done

####Step 1.2: convert to bed file format:
for i in 2 7 ; do
awk -v OFS='\t' '{print $4,$6,$7,$3,$5,$12,$13}' chr${i}.chr2test.blockCoords.txt > chr${i}.chr2test.bed
done

####STEP1.3: Intersect chromosome files
bedtools intersect -wa -wb -a chr2.chr2test.bed -b chr7.chr2test.bed >> intersect.chr2test.bed

```
Testing it on all chromsomes:
```
for{j in Chr01 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09}; do
	grep -E $(echo ${j}) Sb313_vs_ZmB73.blockCoords.txt > ${j}test.blockCoords.txt
	for i in {1..10} ; do
	if grep -Eq $(echo chr${i}) ${j}test.blockCoords.txt
	then
		grep -E $(echo chr${i}) ${j}test.blockCoords.txt > chr${i}.${j}test.blockCoords.txt
		
	fi 
done
done
```

####TRYING TO PARSE RESULTS FROM GENESPACE:

*practicing with sb313 chr10 with B73 and TdFL*

```
#step 1: formatting input files from R to be bedpe (this could be written up in R to generate these files initially)

sed 's/ //g' Sb313_10vsB73.tsv | awk -v OFS='\t' '{print $3,$6,$7,$4,$21,$22,$5,".",".",$20}' - > Sb313_10vsB73.bedpe
sed 's/ //g' Sb313_10vsTdFL.tsv | awk -v OFS='\t' '{print $3,$6,$7,$4,$21,$22,$5,".",".",$20}' - > Sb313_10vsTdFL.bedpe

#step 2: Coverage
#made a bed file for the chromosome 10 in Sb313
tail -n +2 Sb313_10vsB73.bedpe | bedtools coverage -b Sb313_chr10.bed -a -
#It's not doing what I want...
#coverage -d gives me per base coverage...

tail -n +2 Sb313_10vsB73.bedpe | bedtools coverage -d -a Sb313_chr10.bed -b - > Sb313_coverageperbase.txt
awk -v OFS='\t' '$5 != prev; { prev = $5 } ' Sb313_coverageperbase.txt > condensed_Sb313_coverageperbase.txt
#This awk command checks that the previous 5th field != current 5th field; if it doesn't match, then it prints that line
#From this I should be able to make a bed file with the coverages as a range rather than per base
#to make a bed file with these coverage ranges...
tac condensed_Sb313_coverageperbase.txt | awk -v OFS='\t' '{print $1,$4,f,$5} {f=$4-1}' - > Sb313_coverage.bed
#read in the file backwards so that awk can use the starting basepair as the ending basepair from the starting bp of the previous line
awk -v OFS='\t' '$4 ~/1/ || $4~/2/ {print $0}' Sb313_coverage.bed > filtered.Sb313_coverage.bed

#step 3: split the blocks into two temp files that are non-overlapping
tail -n +2 Sb313_10vsB73.bedpe | bedtools intersect -a - -b filtered.Sb313_coverage.bed > Sb313_B73.blocks.intsct.bedpe
cut -f 4 Sb313_B73.blocks.intsct.bedpe | sort | uniq -c | sort -k 1n | tail -n 1 | awk '{print $2}'

for i in chr6 chr5 chr9 ; do grep $i Sb313_B73.blocks.intsct.bedpe > ${i}_Sb313_B73.blocks.intsct.bedpe ; done
#make the intersect files for each B73 chromosome so I can check the overlap
#There should be no overlap between 5 and 9, but there should be overlap with 6 
bedtools window -a chr5_Sb313_B73.blocks.intsct.bedpe -b chr9_Sb313_B73.blocks.intsct.bedpe
bedtools window -a chr5_Sb313_B73.blocks.intsct.bedpe -b chr6_Sb313_B73.blocks.intsct.bedpe

bedtools window -a chr5_Sb313_B73.blocks.intsct.bedpe -b chr6_Sb313_B73.blocks.intsct.bedpe | bedtools overlap -i stdin -cols 2,3,12,13 | awk -v OFS='\t' '$NF !~ /1/ || $NF !~ /-1/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' - >> temp.2
bedtools window -a chr9_Sb313_B73.blocks.intsct.bedpe -b chr6_Sb313_B73.blocks.intsct.bedpe | bedtools overlap -i stdin -cols 2,3,12,13 | awk -v OFS='\t' '$NF !~ /1/ || $NF !~ /-1/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' - >> temp.2

#This gives the chromosome with the maximum number of the blocks

#then the inverse greps
grep -v -f temp.2 chr5_Sb313_B73.blocks.intsct.bedpe >> temp.3
grep -v -f temp.2 chr9_Sb313_B73.blocks.intsct.bedpe >> temp.3
cat temp.3 >> temp.2
rm temp.3

#Then switch the coordinate system to the query system
cut -f 4-10 temp.2 | sort -k1,1 -k2,2n | uniq > temp.2.a
grep -f temp.2.a Sb313_10vsB73.bedpe | sort -k4,4 -k5,5n > temp.2.b

paste temp.2.a temp.2.b | awk -v OFS='\t' '{print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7}' - > Final.B73.2

cut -f 4-10 temp.1 | sort -k1,1 -k2,2n | uniq >temp.1.a
grep -f temp.1.a Sb313_10vsB73.bedpe | sort -k4,4 -k5,5n > temp.1.b

paste temp.1.a temp.1.b |  awk -v OFS='\t' '{print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7}' - > Final.B73.1

#Ok, now let's try to do the TdFL
tail -n +2 Sb313_10vsTdFL.bedpe | bedtools coverage -d -a Sb313_chr10.bed -b - > Sb313_coverageperbase.txt
awk -v OFS='\t' '$5 != prev; { prev = $5 } ' Sb313_coverageperbase.txt > condensed_Sb313_coverageperbase.txt
tac condensed_Sb313_coverageperbase.txt | awk -v OFS='\t' '{print $1,$4,f,$5} {f=$4-1}' - > Sb313_coverage.bed
tail -n +2 Sb313_10vsTdFL.bedpe | bedtools intersect -a - -b filtered.Sb313_coverage.bed > Sb313_TdFL.blocks.intsct.bedpe
cut -f 4 Sb313_TdFL.blocks.intsct.bedpe | sort | uniq -c | sort -k 1n | tail -n 1 | awk '{print $2}'

#chr13, chf4, chr10 (there should be no overlap between 4 and 10)
for i in chr13 chr4 chr10 ; do grep $i Sb313_TdFL.blocks.intsct.bedpe > ${i}_Sb313_TdFL.blocks.intsct.bedpe ; done
bedtools window -a chr4_Sb313_TdFL.blocks.intsct.bedpe -b chr10_Sb313_TdFL.blocks.intsct.bedpe #returns nothing
bedtools window -a chr4_Sb313_TdFL.blocks.intsct.bedpe -b chr13_Sb313_TdFL.blocks.intsct.bedpe 
#put overlaps with chr13 (which should be b) into temp.2
bedtools window -a chr4_Sb313_TdFL.blocks.intsct.bedpe -b chr13_Sb313_TdFL.blocks.intsct.bedpe | bedtools overlap -i stdin -cols 2,3,12,13 | awk -v OFS='\t' '$NF !~ /1/ || $NF !~ /-1/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' - >> temp.2
bedtools window -a chr10_Sb313_TdFL.blocks.intsct.bedpe -b chr13_Sb313_TdFL.blocks.intsct.bedpe | bedtools overlap -i stdin -cols 2,3,12,13 | awk -v OFS='\t' '$NF !~ /1/ || $NF !~ /-1/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' - >> temp.2
cat chr13_Sb313_TdFL.blocks.intsect.bedpe > temp.1

#then the inverse greps
grep -v -f temp.2 chr4_Sb313_TdFL.blocks.intsct.bedpe >> temp.3
grep -v -f temp.2 chr10_Sb313_TdFL.blocks.intsct.bedpe >> temp.3
cat temp.3 >> temp.2
rm temp.3

#Then switch the coordinate system to the query system
cut -f 4-10 temp.2 | sort -k1,1 -k2,2n | uniq > temp.2.a
grep -f temp.2.a Sb313_10vsTdFL.bedpe | sort -k4,4 -k5,5n > temp.2.b

paste temp.2.a temp.2.b | awk -v OFS='\t' '{print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7}' - > temp.2.c

cut -f 4-10 temp.1 | sort -k1,1 -k2,2n | uniq >temp.1.a
grep -f temp.1.a Sb313_10vsTdFL.bedpe | sort -k4,4 -k5,5n > temp.1.b
paste temp.1.a temp.1.b |  awk -v OFS='\t' '{print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7}' - > temp.1.c

#Now different from B73, we have to use the coordinates in B73 vs. TdFL to determin if temp.1 should be final 1 or Final 2
sed 's/ //g' B73vsTdFL.onlySb313Chr10.tsv | awk -v OFS='\t' '{print $3,$6,$7,$4,$21,$22,$5,".",".",$20}' - > B73vsTdFL_Sb313.10.bedpe
#NOTE that we'll need a way to know if B73 or query is chr1 or chr2 when automating
#manually we know TdFL is genome 1 and B73 is genome 2
#issue with the temp.2.c and intersect because of the "-" and the end being earlier than the start bp
awk -v OFS='\t' '{ if ($3 < $2 ) print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10 ; else print $0 } ' temp.2.c > temp.2.d
awk -v OFS='\t' '{ if ($3 < $2 ) print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10 ; else print $0 } ' temp.1.c > temp.1.d

tail -n +2 B73vsTdFL_Sb313.10.bedpe | bedtools intersect -a temp.2.d -b -
awk -v OFS='\t' '{if ($3 < $2) print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10; else print $0 }' Final.B73.1 > edited.Final.B73.1
awk -v OFS='\t' '{if ($3 < $2) print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10; else print $0 }' Final.B73.2 > edited.Final.B73.2

#This is to find which B73 chr our temp2 sequences for TdFL match to
#note there is one match from TdFLchr4 to ZmB73chr6 while all other TdFLchr4 match to B73chr9
tail -n +2 B73vsTdFL_Sb313.10.bedpe | bedtools intersect -a - -b temp.2.d -wb

#Since we expect that no single chromosome should have both subgenome 1 and subgenome 2 on it... 
#but we won't know which is which...
tail -n +2 B73vsTdFL_Sb313.10.bedpe | bedtools intersect -a - -b temp.2.d -wb > temp.2.check.bedpe
tail -n +2 B73vsTdFL_Sb313.10.bedpe | bedtools intersect -a - -b temp.1.d -wb > temp.1.check.bedpe
bedtools intersect -a temp.1.check.bedpe -b temp.2.check.bedpe 
#This shows that there's no overlapping of the TdFL coordinates (which is good)
 
tail -n +2 B73vsTdFL_Sb313.10.bedpe | bedtools intersect -a - -b temp.2.d -wb | awk -v OFS='\t' '{if($6 > $5) print $4,$5,$6,$1,$2,$3 ; else print $4,$6,$5,$1,$2,$3}' - > temp.2.check2.bedpe
tail -n +2 B73vsTdFL_Sb313.10.bedpe | bedtools intersect -a - -b temp.1.d -wb | awk -v OFS='\t' '{if($6 > $5) print $4,$5,$6,$1,$2,$3 ; else print $4,$6,$5,$1,$2,$3}' - > temp.1.check2.bedpe

bedtools intersect -a temp.1.check2.bedpe -b temp.2.check2.bedpe -wa -wb
#this shows that there are 2 TdFL coordinates for one spot of B73
#Do we need to remove one set? Do we need just assign them to the file with the most like? 
#Do I go by which file has the greatest number of bps matching to assign final 1 and 2?
#I think we want to keep them both. We'll just assign the problematic (?) one to the file 1 or 2 where the rest of that chromosome has been assigned?
 
 
bedtools intersect -a edited.Final.B73.2 -b temp.2.check2.bedpe -wb |cut -f 11-16 | sort | uniq | awk -v OFS='\t' '$4 ~ /chr4/ {print $0}' - > temp.2.vsB73.2.bedpe
bedtools intersect -a edited.Final.B73.1 -b temp.2.check2.bedpe -wb |cut -f11-16 | sort | uniq | awk -v OFS='\t' '$4 ~ /chr4/ {print $0}' -> temp.2.vsB73.1.bedpe

bedtools intersect -a edited.Final.B73.1 -b temp.2.check2.bedpe -wb |cut -f11-16 | sort | uniq | awk -v OFS='\t' '$4 ~ /chr4/ {print $6-$5}' - > temp.2.vsB73.1.sum

awk '{sum += $1 } END {print sum}' temp.2.vsB73.2.sum
```
Let's try to automate the above:

```
#!/bin/bash

QUERY=$1 #B73
REFCHR=$2 #10

module load bedtools2

#step 1: formatting input files from R to be bedpe (this could be written up in R to generate these files initially)

sed 's/ //g' Sb313_${REFCHR}vs${QUERY}.tsv | awk -v OFS='\t' '{print $3,$6,$7,$4,$21,$22,$5,".",".",$20}' - | tail -n +2 > Sb313_${REFCHR}vs${QUERY}.bedpe

#step 2: Coverage 
####NOTE: NEED TO MAKE CHR BED FILES FOR EACH SB313 FROM THE LENGTH IN THE FAI; NEED A NEW NAME FOR WRITTEN FILES FROM THIS POINT ONWARD IF TO BE RUN PARALLEL

bedtools coverage -d -a Sb313_chr${REFCHR}.bed -b Sb313_${REFCHR}vs${QUERY}.bedpe > Sb313_coverageperbase.txt
awk -v OFS='\t' '$5 != prev; { prev = $5 } ' Sb313_coverageperbase.txt > condensed_Sb313_coverageperbase.txt
#This awk command checks that the previous 5th field != current 5th field; if it doesn't match, then it prints that line
#From this I should be able to make a bed file with the coverages as a range rather than per base

#to make a bed file with these coverage ranges...
tac condensed_Sb313_coverageperbase.txt | awk -v OFS='\t' '{print $1,$4,f,$5} {f=$4-1}' - > Sb313_coverage.bed
#read in the file backwards so that awk can use the starting basepair as the ending basepair from the starting bp of the previous line

#only keep regions of Sb313 REFCHR that have coverage == 1 or == 2
awk -v OFS='\t' '$4 ~/1/ || $4~/2/ {print $0}' Sb313_coverage.bed > filtered.Sb313_coverage.bed

#step 3: split the blocks into two temp files that are non-overlapping

#Finds intersection of the block coordinates and the regions of REFCHR that only have coverage == 1 or ==2
bedtools intersect -a Sb313_${REFCHR}vs${QUERY}.bedpe -b filtered.Sb313_coverage.bed > Sb313_${QUERY}.blocks.intsct.bedpe
#Find the QUERY CHR that has the most blocks from the intersection
#This chromosome will be written to temp.1
cut -f 4 Sb313_${QUERY}.blocks.intsct.bedpe | sort | uniq -c | sort -k 1n | tail -n 1 | awk '{print $2}'

#NOTE: I need to get from the previous arguments to a list of the QUERY CHRs 
#NOTE: THIS LIST NEEDS TO BE ABLE TO BE USED LIKE ARGUMENTS IN THE FUTURE STEPS
#Here that is chr6 (which is the max number of blocks that should go to temp.1), chr 5, and chr9

for i in chr6 chr5 chr9 ; do grep $i Sb313_${QUERY}.blocks.intsct.bedpe > ${i}_Sb313_${QUERY}.blocks.intsct.bedpe ; done

cat chr6_Sb313_${QUERY}.blocks.intsect.bedpe > temp.1

#Find the intersection of the remaining query chromosomes with the one that has already been assigned to temp.1
bedtools window -a chr5_Sb313_B73.blocks.intsct.bedpe -b chr6_Sb313_B73.blocks.intsct.bedpe | bedtools overlap -i stdin -cols 2,3,12,13 | awk -v OFS='\t' '$NF !~ /1/ || $NF !~ /-1/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' - >> temp.2
bedtools window -a chr9_Sb313_B73.blocks.intsct.bedpe -b chr6_Sb313_B73.blocks.intsct.bedpe | bedtools overlap -i stdin -cols 2,3,12,13 | awk -v OFS='\t' '$NF !~ /1/ || $NF !~ /-1/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' - >> temp.2

#then the inverse greps to get the remaining blocks of those chromosomes that don't overlap with chr6
#ASSUMES THAT THEY DON'T OVERLAP WITH EACH OTHER AT ALL; WHICH IS A FINE ASSUMPTION WHEN THERE ARE 3 QUERY CHROMOSOMES, BUT BREAKS WHEN THERE'S MORE THAN 3
grep -v -f temp.2 chr5_Sb313_B73.blocks.intsct.bedpe >> temp.3
grep -v -f temp.2 chr9_Sb313_B73.blocks.intsct.bedpe >> temp.3
cat temp.3 >> temp.2
rm temp.3

#NOTE: THIS SHOULD ONLY BE IF THE QUERY == B73
#Then switch the coordinate system to the query system
cut -f 4-10 temp.2 | sort -k1,1 -k2,2n | uniq > temp.2.a
grep -f temp.2.a Sb313_${REFCHR}vs${QUERY}.bedpe | sort -k4,4 -k5,5n > temp.2.b
#if statement reverses coordinates if the end < start bp, which happens if orientation (strand) is -
paste temp.2.a temp.2.b | awk -v OFS='\t' '{if ($3 > $2) print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7 ; else print $1,$3,$2,$8,$9,$10,$4,$5,$6,$7}' - > Final.${QUERY}.2

#Then do the same for temp.1
cut -f 4-10 temp.1 | sort -k1,1 -k2,2n | uniq >temp.1.a
grep -f temp.1.a Sb313_10vsB73.bedpe | sort -k4,4 -k5,5n > temp.1.b
paste temp.1.a temp.1.b |  awk -v OFS='\t' '{if ($3 > $2) print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7 ; else print $1,$3,$2,$8,$9,$10,$4,$5,$6,$7}' - > Final.${QUERY}.1

#NOTE: IF THE QUERY != B73 THEN WE HAVE TO DO AN EXTRA CHECK TO SEE IF IT'S TEMP1 SHOULD BE FINAL 1 OR FINAL 2

#NOTE: WHILE BY THIS POINT ALL QUERY TO SB313 SHOULD BE THE APPROPRIATE COVERAGE (1 OR 2)
#BUT THERE MAY BE POINTS WHERE THE QUERY IS 2:1 B73; WE DON'T CARE IF IT'S NOT 1:1
#BUT THAT DOES MUDDY THE WATERS WHEN DECIDING IF TEMP1 SHOULD BE FINAL 1 OR FINAL 2
#SO IF THAT HAPPENS WHERE 1 QUERY COORD. MATCHES B73 FINAL FILE 1 AND THE OTHER MATCHES B73 FINAL FILE 2
#THE QUERY CHR SHOULD GO TO THE FILE WHERE THE OTHER BLOCKS ON THAT QUERY CHROMOSOME MATCH THE MOST
#EX: CHR 4 TDFL MATCHES MOSTLY TO CHR 9 IN B73 BUT THERE'S A ~2 MILLION BP CHUNK THAT MATCHES TO CHR 6 IN B73
#THIS SHOULD BE RESOLVED WITH CHR 4 BLOCKS GOING TO FINAL FILE 2 (MATCHING B73 CHR9'S PLACEMENT) INSTEAD OF SPLITTING CHR4 ACROSS FINAL FILES 1 AND 2

#Now different from B73, we have to use the coordinates in B73 vs. TdFL to determine if temp.1 should be final 1 or Final 2
sed 's/ //g' B73vs${QUERY}.onlySb313Chr${REFCHR}.tsv | awk -v OFS='\t' '{print $3,$6,$7,$4,$21,$22,$5,".",".",$20}' -|tail -n +2 > B73vs${QUERY}_Sb313.${REFCHR}.bedpe
#NOTE: WE'LL NEED A CHECK TO SEE IF B73 IS GENOME 1 OR GENOME 2
#manually we know TdFL is genome 1 and B73 is genome 2

cut -f 4-10 temp.2 | sort -k1,1 -k2,2n | uniq > temp.2.a
grep -f temp.2.a Sb313_${REFCHR}vs${QUERY}.bedpe | sort -k4,4 -k5,5n > temp.2.b
#if statement reverses coordinates if the end < start bp, which happens if orientation (strand) is -
paste temp.2.a temp.2.b | awk -v OFS='\t' '{if ($3 > $2) print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7 ; else print $1,$3,$2,$8,$9,$10,$4,$5,$6,$7}' - > temp.2.c

#Then do the same for temp.1
cut -f 4-10 temp.1 | sort -k1,1 -k2,2n | uniq >temp.1.a
grep -f temp.1.a Sb313_10vsB73.bedpe | sort -k4,4 -k5,5n > temp.1.b
paste temp.1.a temp.1.b |  awk -v OFS='\t' '{if ($3 > $2) print $1,$2,$3,$8,$9,$10,$4,$5,$6,$7 ; else print $1,$3,$2,$8,$9,$10,$4,$5,$6,$7}' - > temp.2.c

#Now check to see if temp.c.1 or temp.c.2 are best matched to Final.B73.1 or Final.B73.2

bedtools intersect -a B73vs${QUERY}_Sb313.${REFCHR}.bedpe -b temp.2.c -wb > temp.2.check.bedpe
bedtools intersect -a B73vs${QUERY}_Sb313.${REFCHR}.bedpe -b temp.1.c -wb > temp.2.check.bedpe
bedtools intersect -a temp.1.check.bedpe -b temp.2.check.bedpe 
#This should show that there's no overlapping of the QUERY coordinates between temp files 1 and 2 (which is good)

#switch to B73 COORDINATES VS. QUERY
awk -v OFS='\t' '{if($6 > $5) print $4,$5,$6,$1,$2,$3 ; else print $4,$6,$5,$1,$2,$3}' temp.2.check.bedpe > temp.2.check2.bedpe
awk -v OFS='\t' '{if($6 > $5) print $4,$5,$6,$1,$2,$3 ; else print $4,$6,$5,$1,$2,$3}' temp.1.check.bedpe > temp.1.check2.bedpe

#This shows the 1 overlap in the test of B73 and TdFL
bedtools intersect -a Final.B73.1 -b temp.2.check.bedpe
bedtools intersect -a Final.B73.2 -b temp.2.check.bedpe

###NOTE: HOW TO TELL IT BASED OFF THESE INTERSECT RUNS WHICH SHOULD GO WHICH WITHOUT MANUAL INSPECTION?
#FROM MANUAL INSPECTION FOR B73 VS TDFL

cat temp.2.c > Final.${QUERY}.2
cat temp.1.c > Final.${QUERY}.1
```
Attempting Arun's script
```
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:nova21-38:33421 snodgras@nova.its.iastate.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: snodgras
   password: JW5E+Rs3XzTJGDbvDEa1

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f 4687687
```

*ARUN RAN HIS SCRIPT AND CREATED THE BINNING*
*I'M MOVING THE RESULTS FROM MY LOCAL COMPUTER TO NOVA*

#Run AnchorWave on Complete Genomes
##Minimap
Script: `slurm_02.minimapAnchorPoints.sh` and `02.minimapAnchorPoints.sh` (the if statement has weird colors in vi, not sure if it'll run...)
Time to Run: 00:01:29:00 ish

Make a directory for the sam files to go to: 
```
mkdir sam_files
```
And make a softlink for the Av genome to the split genomes file
```
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta split_genome_assemblies/Av.fasta
```
```bash script
FILE1=$1 #Genomic fasta file ID for maize genome (genomeID)) (example: Zm-B73-REFERENCE-NAM-5.0); be sure to make certain the fasta file extension in the script matches user file extension i.e. '.fasta' vs '.fa', etc)
FileName=${FILE1#split_genome_assemblies/}
ml minimap2

#Do I need to redo the minimap to sorghum to itself everytime? Could cut down with a if file exists exemption
if [[ ! -f "sam_files/Sb313.ref.sam" ]] ; then
	minimap2 \
		-x splice \
		-t 11 \
		-k 12 \
		-a \
		-p 0.4 \
		-N 20 \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > sam_files/Sb313.ref.sam
fi

minimap2 \
	-x splice \ #Applies multiple options at once: long-read spliced alignment; long deletions are taken as introns and represented as "n" in CIGAR, long insertions are disabled, indel gap costs are different during chaining; computing 'ms' tag ignores introns to demote hits to pseudogenes
	-t 11 \ #threads
	-k 12 \ #minimizer k-mer length
	-a \ #generate CIGAR and ouput alignmnets in SAM format
	-p 0.4 \ #Minimal secondary-to-primary score ratio to output secondary mappings
	-N 20 \ #Output at most INT secondary alignments
	${FILE1} \ #target fasta
	Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > sam_files/${FileName%.fasta}_Sb313.cds.sam #query fasta and output file
```

```slurm script
for i in split_genome_assemblies/*.fasta ; do
        bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/02.minimapAnchorPoints.sh ${i}
;done
```
##Run Anchorwave
Make sure that Av is 1:1 and the rest is 2:1
```03.runAnchorWave.Tripsacinae.sh
#!/bin/bash
FILE1=$1
FileName=${FILE1#*/}
GenomeName=$2
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \ 
	-as Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \ 
	-a sam_files/${FileName}_Sb313.cds.sam \ 
	-ar sam_files/Sb313.ref.sam \ 
	-s ${FILE1}.fasta \ 
	-n AnchorWave_output/Sb313_${GenomeName}_anchorwave.anchors \ 
	-R 2 \ 
	-Q 1 \
	-o AnchorWave_output/Sb313_${GenomeName}_anchorwave.maf \ 
	-t 5 > AnchorWave_logs/Sb313_${GenomeName}_anchorwave.log 2>&1 
	
03.runAnchorWave.Av.sh
#!/bin/bash

FILE1=$1
FileName=${FILE1#*/}
GenomeName=$2
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \
	-as Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
	-a sam_files/${FileName}_Sb313.cds.sam \
	-ar sam_files/Sb313.ref.sam \
	-s ${FILE1}.fasta \
	-n AnchorWave_output/Sb313_${FileName}_anchorwave.anchors \
	-R 1 \
	-Q 1 \
	-o AnchorWave_output/Sb313_${FileName}_anchorwave.maf \
	-t 5 > AnchorWave_logs/Sb313_${FileName}_anchorwave.log 2>&1

```
```
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a TdFL >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Av.sh assemblies_final/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0 Av  >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zn-PI615697-REFERENCE-PanAnd-1.0 ZnPI615697 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0 ZvTIL01 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zd-Gigi-REFERENCE-PanAnd-1.0 ZdGigi >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zv-TIL11-REFERENCE-PanAnd-1.0 ZvTIL11 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zd-Momo-REFERENCE-PanAnd-1.0 ZdMomo >> scripts/AnchorWave.cmds.txt    
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0 ZxTIL18 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zh-RIMHU001-REFERENCE-PanAnd-1.0 ZhRIMHU001 >> scripts/AnchorWave.cmds.txt      
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zx-TIL25-REFERENCE-PanAnd-1.0 ZxTIL25 >> scripts/AnchorWave.cmds.txt



echo bash scripts/03.runAnchorWave.Tripsacinae.sh
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-B73 ZmB73
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-B97 ZmB97
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML103 ZmCML103
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML228 ZmCML228
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML247 ZmCML247
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML277 ZmCML277
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML322 ZmCML322
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML333 ZmCML333
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML52 ZmCML52
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML69 ZmCML69
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-HP301 ZmHP301
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-IL14H ZmIL14H
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Ki11 ZmKi11
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Ki3 ZmKi3
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Ky21 ZmKy21
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-M162W ZmM162W
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-M37W ZmM37W
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Mo18W ZmMo18W
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-MS71 ZmMS71
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-NC350 ZmNC350
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-NC358 ZmNC358
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Oh43 ZmOh43
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Oh7b ZmOh7b
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-P39 ZmP39
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Tx303 ZmTx303
bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Tzi8 ZmTzi8 >> scripts/AnchorWave.cmds.txt

python makeSLURM.py 1 AnchorWave.cmds.txt

```
For testing, to cut down on time, we're testing this out on a subset of lines
```
for i in {0..4} ; do sbatch scripts/AnchorWave.cmds_${i}.sub ; done
for i in {10..19} ; do sbatch scripts/AnchorWave.cmds_${i}.sub ; done
```
After testing was complete, I ran the other AnchorWave commands:
```
for i in {5..9}; do sbatch scripts/AnchorWave.cmds_${i}.sub ; done
for i in {20..35} ; do sbatch scripts/AnchorWave.cmds_${i}.sub ; done
```
Want to include luxurians both in alignment to Sorghum and to B73
```
#assemblies-final/Zl-RIL003-REFERENCE-PanAnd-1.0.fasta
cp AnchorWave.cmds_1.sub AnchorWave.cmds_36.sub
#edit to be 
bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zl-RIL003-REFERENCE-PanAnd-1.0 ZlRIL003

#then make separate slurm script that will use B73 as reference
#slurm_03.ZmB73vsZlux.AnchorWave.sh
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
	assemblies_final/Zl-RIL003-REFERENCE-PanAnd-1.0 Zl-RIL003.fasta  \
	NAM-assemblies/Zm-B73.fasta  > sam_files/ZlRIL003_Sb313.cds.sam

ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i NAM-annotations/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 \
	-as NAM-assemblies/Zm-B73.fasta \
	-r NAM-assemblies/Zm-B73 \
	-a sam_files/ZlRIL003_ZmB73.cds.sam \
	-ar sam_files/ZmB73.ref.sam \
	-s assemblies_final/Zl-RIL003-REFERENCE-PanAnd-1.0 Zl-RIL003.fasta \
	-n AnchorWave_output/ZmB73_ZlRIL003_anchorwave.anchors \
	-R 1 \
	-Q 1 \
	-o AnchorWave_output/ZmB73_ZlRIL003_anchorwave.maf \
	-t 5 > AnchorWave_logs/ZmB73_ZlRIL003_anchorwave.log 2>&1 

```

#Swap coordinates of MAF
Use script `https://github.com/baoxingsong/AnchorWave/blob/master/scripts/anchorwave-maf-swap.py`
```
genome=$1
cat Sb313_${genome}_anchorwave.maf | python /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/anchorwave-maf-swap.py > Sb313_${genome}_anchorwave.swap.maf

###SLURM
bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/03.5.swapMAF.sh TdFL
```

#Split MAF by parsing coordinates
_Check MAF with IGV_
I can't get IGV to load MAF Correctly: `https://github.com/igvteam/igv/blob/master/test/data/maf_ucsc/ucscSample.maf`
*Tool Options*
`/ptmp/arnstrm/maftools.sif`
mafExtractor 
A program to extract all alignment blocks that contain a region in a particular sequence. 
Useful for isolating regions of interest in large maf files.
`https://github.com/dentearl/mafTools`
It does change the `a=score` to 0 and adds `mafExtractor_splicedBlock=true splice_id=1_0`

```
while read line ; do singularity exec maftools.sif mafExtractor --maf AnchorWave_output/Sb313_TdFL_anchorwave.maf --seq $(echo $line | cut -f 1 -d " ") --start $(echo $line | cut -f 2 -d " ") --stop $(echo $line |cut -f 3 -d " ") --soft >> test.maf ; done < test.txt 
```
*This gives me two regions instead of 1*
```
ml singularity
singularity exec --bind $PWD maftools.sif mafExtractor\
	--maf maf \
	--seq \
	--start \
	--stop \
	 > split.maf
```

`/work/LAS/mhufford-lab/maffilter`
`https://jydu.github.io/maffilter/Manual/SelectChr.html#SelectChr`
maffilter
*This has no way of writing the filtered output as a file? Where does the output even goooo???*
```
#Change the seqnames to include the genome names
sed -i "s/Chr/Sb313.Chr/g" Sb313_TdFL_tabs.maf
sed -i "s/chr/TdFL.chr/g" Sb313_TdFL_tabs.maf
sed -i "s/scaf_B/TdFL.scaf_B/g" Sb313_TdFL_tabs.maf

/work/LAS/mhufford-lab/maffilter param=maffilter.optionfile

maf.filter=SelectChr(ref_species=species1, chromosome=chr1)
ExtractFeature(ref_species=TdFL, feature.file=parsing_coordinates/uniq.zmB73vsTdFL.bin1.bed, feature.format=BedGraph)
#Also above doesn't work because bedGraph != bed
```

*Manual Options*

```
MAF=$1
genome=$2

grep -w "chr" $MAF | cut -f 2-5 -d " " > temp.coords

#awk -v OFS='\t' '$4 ~ /+/ {print $1, $2, $2+$3}' temp.coords >> temp.bed #positive strand
#awk -v OFS='\t' '$4 ~ /-/ {print $1, $2-$3, $2}' temp.coords >> temp.bed #negative strand

awk -v OFS='\t' '{print $1, $2, $2+$3, $4}' temp.coords >> temp.bed #incase strand is already accounted for in MAF coordinates

bedtools intersect temp.bed bin1.parsing.bed > temp.bin1.bed
bedtools intersect temp.bed bin2.parsing.bed > temp.bin2.bed

awk -v OFS='\t' '{print "s", $1,$2,$3-$2,$4}' temp.bin1.bed > temp.bin1.coords
awk -v OFS='\t' '{print "s", $1,$2,$3-$2,$4}' temp.bin2.bed > temp.bin2.coords

grep -B2 -A1 -w -f temp.bin1.coords $MAF >> ${genome}.bin1.MAF
grep -B2 -A1 -w -f temp.bin2.coords $MAF >> ${genome}.bin2.MAF
```
Testing:
```
grep "chr" AnchorWave_output/Sb313_TdFL_anchorwave.maf | cut -f 2-5 > temp.coords

awk '$4 ~ /-/ {print $2-$3}' temp.coords #negative numbers indicate that maf has already taken strandedness into account?

awk -v OFS='\t' '{print $1, $2, $2+$3, $4}' temp.coords >> temp.bed #incase strand is already accounted for in MAF coordinates

ml bedtools2
bedtools intersect -a temp.bed -b parsing_coordinates/uniq.zmB73vsTdFL.bin1.bed > temp.bin1.bed

#looks like there's duplicates from the intersect:
grep "chr2    139405908       140111349" temp.bed 
grep "chr2    139405908       140111349" temp.bin1.bed 

sort temp.bin1.bed | uniq > uniq.temp.bin1.bed

bedtools intersect -a temp.bed -b parsing_coordinates/uniq.zmB73vsTdFL.bin2.bed > temp.bin2.bed
sort temp.bin2.bed | uniq > uniq.temp.bin2.bed

#check that they're exclusive:
bedtools intersect -wa -wb -a uniq.temp.bin1.bed -b uniq.temp.bin2.bed

bedtools intersect -wa -a uniq.temp.bin1.bed -b uniq.temp.bin2.bed > bin1.dups.bed
bedtools intersect -wb -a uniq.temp.bin1.bed -b uniq.temp.bin2.bed > bin2.dups.bed

grep -v -f bin1.dups.bed uniq.temp.bin1.bed > excl.uniq.bin1.bed
grep -v -f bin2.dups.bed uniq.temp.bin2.bed > excl.uniq.bin2.bed

awk -v OFS='\t' '{print "s",$1,$2,$3-$2,$4}' excl.uniq.bin1.bed > temp.bin1.coords
awk -v OFS='\t' '{print "s",$1,$2,$3-$2,$4}' excl.uniq.bin2.bed > temp.bin2.coords

awk -v OFS='\t' '{print $0}' AnchorWave_output/Sb313_TdFL_anchorwave.maf > Sb313_TdFL_tabs.maf

#for manual inspection: awk -v OFS='\t' -v c="chr10" -v s="17816503" "'$2 == c && $3 == s {print g; print f; print $1,$2,$3,$4,$5; printf "\n"} {g=f} {f=$2}' Sb313_TdFL_tabs.maf 
#awk -v OFS='\t' '$2 == "chr10" && $3 == 17816503 {print g; print f; print $1,$2,$3,$4,$5; printf "\n"} {g=f} {f=$2}' AnchorWave_output/Sb313_TdFL_anchorwave.maf

awk -v OFS='\t' '$2 == "chr10" && $3 == 17816503 {print g; print f; print $0; printf "\n"} {g=f} {f=$0}' Sb313_TdFL_tabs.maf 
```
*Can't loop through Awk*
*NOTE: There is no chr6 in the bin coordinates for TdFL*

*kentutils*
`/ptmp/arnstrm/kentutils.sif`

see `split_and_filter_maf_files_for_Sam.txt`

```mafsplit.1.sh
#!/bin/bash
ml singularity
genome=$1
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName swap_${genome}/ Sb313_${genome}_anchorwave.swap.maf

###
for i in TdFL ZdGigi ZvTIL01 ZmB73 ; do echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/mafsplit.1.sh $i >> ../scripts/mafsplit.cmds.txt ; done
cd ../scripts
python makeSLURM.py 1 mafsplit.cmds.txt
cd -
sbatch ../scripts/mafsplit.cmds_0.sub
for i in {1..3} ; do sbatch --dependency=afterok:4749635 ../scripts/mafsplit.cmds_${i}.sub ; done

###

for f in */*.maf ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp"_"$fl"; done
rm *scaf*
```

```mafsplit.2.sh
#!/bin/bash

ml python/2.7.18-2ut3ogj #or however python2 is invoked in Nova, if necessary
#ml kentutils #or whichever module has this tool in Nova
ml singularity

for input in swap_*.maf
  do
        echo $input
        output=$(echo ${input} | sed 's/.maf//g; s/swap_//g')
        echo $output

cat ${input} | python /scripts/anchorwave-maf-swap.py  > ${output}_split.maf


#STEP 4:
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif  mafSplit -byTarget dummy.bed -useFullSequenceName ${output}/ ${output}_split.maf

done

for f in */*.maf ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp"_"$fl"; 

done
```

Manually:
```
mkdir zm-sb_split_mafs
mv *chr*_Chr*.maf zm-sb_split_mafs
cd zm-sb_split_mafs
```

Next, sort maf files, then remove overlapping blocks that cause Zack's gvcf converter to fail 
(these aren't homeolog repeats since everything is split by homeolog):
maf-sort.sh: get the code for this from here: 
https://github.com/UCSantaCruzComputationalGenomicsLab/last/blob/master/scripts/maf-sort.sh

```mafsplit.3.sh 

#!/bin/bash

ml singularity

for input in *.maf
  do
      echo $input
        output=$(echo ${input} | sed 's/.maf//g')
        echo $output

sh /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/maf-sort.sh  ${input} > ${output}_sorted.maf

singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafFilter ${output}_sorted.maf -overlap > ${output}_sorted_filtered.maf

done
```
merge split mafs into bins for chr10 test:
```
cat TdFL_chr13_Chr10_sorted_filtered.maf > TdFL_Chr10.bin2.maf
cat TdFL_chr10_Chr10_sorted_filtered.maf > TdFL_Chr10.bin1.maf
tail -n +2 TdFL_chr4_Chr10_sorted_filtered.maf >> TdFL_Chr10.bin1.maf

cat ZdGigi_chr5_Chr10_sorted_filtered.maf  > ZdGigi_Chr10.bin1.maf
tail -n +2 ZdGigi_chr9_Chr10_sorted_filtered.maf  >> ZdGigi_Chr10.bin1.maf
cat ZdGigi_chr6_Chr10_sorted_filtered.maf  > ZdGigi_Chr10.bin2.maf

cat ZvTIL01_chr5_Chr10_sorted_filtered.maf  > ZvTIL01_Chr10.bin1.maf
tail -n +2 ZvTIL01_chr9_Chr10_sorted_filtered.maf  >> ZvTIL01_Chr10.bin1.maf
cat ZvTIL01_chr6_Chr10_sorted_filtered.maf  > ZvTIL01_Chr10.bin2.maf

cat ZmB73_chr5_Chr10_sorted_filtered.maf  > ZmB73_Chr10.bin1.maf
tail -n +2 ZmB73_chr9_Chr10_sorted_filtered.maf  >> ZmB73_Chr10.bin1.maf
cat ZmB73_chr6_Chr10_sorted_filtered.maf  > ZmB73_Chr10.bin2.maf

wc -l *Chr10_sorted_filtered.maf #668
wc -l *Chr10.bin*maf #664 (4 less because removing header line of 1 maf file)
```

Try converting to GVCF and see if it works
```
for i in *bin*.maf ; do 
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/04.splitmaf2gvcf.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.fasta ${i} >> ../../scripts/splitmaf2gvcf.cmds.txt
done
cd ../../scripts/
python makeSLURM.py 1 splitmaf2gvcf.cmds.txt
cd - 
#test before running all
sbatch ../../scripts/splitmaf2gvcf.cmds_0.sub 

for i in {1..7} ; do sbatch ../../scripts/splitmaf2gvcf.cmds_${i}.sub ; done
```
Trimming the GVCFs
```
for i in *gvcf.gz ; do 
echo module load samtools \; module load gatk \; gatk LeftAlignAndTrimVariants -O trimmed_${i} -V ${i} -R /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.fasta --max-indel-length 9101263 &> gf_stdout.txt >> ../../scripts/splitmaf.trimGVCF.cmds.txt
done
cd -
#edit makeSLURM.py to be only 1 hour of wall time
python makeSLURM.py 1 splitmaf.trimGVCF.cmds.txt
cd -
sbatch ../../scripts/splitmaf.trimGVCF.cmds_0.sub
for i in {1..7} ; do sbatch --dependency=afterok: ../../scripts/splitmaf.trimGVCF.cmds_${i}.sub ;done
```
Then manually remove the too large indels
```slurm_04.5.removelargeindels.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=2:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml vcftools
awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_0.e4756881 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_TdFL_Chr10.bin1.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_TdFL_Chr10.bin1.gvcf.gz

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_1.e4756882 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_TdFL_Chr10.bin2.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_TdFL_Chr10.bin2.gvcf.gz

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_2.e4756883 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZdGigi_Chr10.bin1.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZdGigi_Chr10.bin1.gvcf.gz

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_3.e4756884 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZdGigi_Chr10.bin2.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZdGigi_Chr10.bin2.gvcf.gz

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_4.e4756885 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZmB73_Chr10.bin1.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZmB73_Chr10.bin1.gvcf.gz

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_5.e4756886 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZmB73_Chr10.bin2.gvcf.gz  --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZmB73_Chr10.bin2.gvcf.gz 

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_6.e4756887 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZvTIL01_Chr10.bin1.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZvTIL01_Chr10.bin1.gvcf.gz

awk '$0 ~ /Indel is too long/  { print }' splitmaf.trimGVCF.cmds_7.e4756888 |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZvTIL01_Chr10.bin2.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZvTIL01_Chr10.bin2.gvcf.gz 

ml gatk
for i in ready*recode.vcf ; do gatk IndexFeatureFile -F $i ;done

```
output is like: `ready_Sb313_TdFL_Chr10.bin1.gvcf.gz.recode.vcf`

convert to vcf
```slurm_06.10.splitmaf.gvcf2vcf.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=96:00:00
#SBATCH --job-name=gatk-chr10
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml gatk
ref=/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.chr10.fasta
ml samtools
# index
samtools faidx $ref
# dict
gatk CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${ref%.*}.dict
# db import
gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
-V ready_Sb313_TdFL_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_TdFL_Chr10.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdGigi_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdGigi_Chr10.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB73_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB73_Chr10.bin2.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL01_Chr10.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL01_Chr10.bin2.gvcf.gz.recode.vcf \
--batch-size 1 \
--genomicsdb-workspace-path chr1_gatkDBimport \
-L 10 \
--genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 
#--tmp-dir $TMPDIR

# vcf output
gatk --java-options "-Xmx50g" GenotypeGVCFs \
-R $ref \
-V gendb://chr10_gatkDBimport \
--cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000 \
-O chr10.vcf
```
prep reference annotation to be included as a filter for vcf
```slurm_07.00.splitmaf.genebed.sh
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
bedops gff2bed < /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 > /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.gene.bed
```
^^^ Didn't work so I'll try manually making a bed file
```
#note there are no exon features in the Sb313 gff3

awk -v OFS='\t' '$3 == "CDS" && $7 == "+" {print $1,$4-1,$5,$9,$7}; $3 == "CDS" && $7 == "-" {print $1,$4,$5+1,$9,$7}' Sbicolor_313_v3.1.gene.gff3 > Sb313.CDS.bed
sed -i "s/Chr//g" Sb313.CDS.bed

awk -v OFS='\t' '$3 == "gene" && $7 == "+" {print $1,$4-1,$5,$9,$7}; $3 == "CDS" && $7 == "-" {print $1,$4,$5+1,$9,$7}' Sbicolor_313_v3.1.gene.gff3 > Sb313.gene.bed
sed -i "s/Chr//g" Sb313.gene.bed
```
then subset out the indels
```slurm_07.10.splitmaf.subsetvcf.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=4:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
ml gatk
#gatk SelectVariants -V chr10.vcf -O chr10.indelonly.vcf --select-type-to-include INDEL
gatk SelectVariants -V chr10.vcf -O chr10.CDSonly.indelonly.vcf --select-type-to-include INDEL -L /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb313.CDS.bed 
gatk SelectVariants -V chr10.vcf -O chr10.geneonly.indelonly.vcf --select-type-to-include INDEL -L /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb313.gene.bed
```
then convert to a table? Or some other way of intersecting with Sb exon coordinates? 
```slurm_08.10.splitmaf.vcftable.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=4:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#ml gatk
#gatk VariantsToTable -V chr10.indelonly.vcf -O chr10.indelonly.table -F CHROM -F POS -F TYPE -F NCALLED -GF GT -GF GQ
#gatk VariantsToTable -V chr10.CDSonly.indelonly.vcf -O chr10.CDS.only.indelonly.table -F CHROM -F POS -F TYPE -F NCALLED -GF GT -GF GQ
ml bcftools
bcftools annotate -x INFO,^FORMAT/GT chr10.CDSonly.indelonly.vcf  > chr10.CDSonly.indelonly_reformted.vcf
bcftools annotate -x INFO,^FORMAT/GT chr10.geneonly.indelonly.vcf  > chr10.geneonly.indelonly_reformatted.vcf

```
Make a key for parsing in R
```
awk -v OFS='\t' '$1 == 10 {print $1,$2,$4,$5}' chr10.CDSonly.indelonly.vcf > chr10.CDSonly.indelonly.refaltKey.tsv
```
Use `SNPeff` to add annotations
```
java -Xmx8g -jar /work/LAS/mhufford-lab/snodgras/Fractionation/snpEff/snpEff.jar \
-interval /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb313.CDS.bed \
chr10.CDSonly.indelonly_reformted.vcf > chr10.CDSonly.indelonly.refmt.ann.vcf

ml snpeff
snpEff -interval /work/LAS/mhufford-lab/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sb313.CDS.bed \
chr10.CDSonly.indelonly_reformted.vcf > chr10.CDSonly.indelonly.refmt.ann.vcf
```

For using Rstudio on nova:
```
grep -v "^##" chr10.geneonly.indelonly_reformatted.vcf > chr10.geneonly.indelonly_reformatted.tsv
 ssh -N -L 8787:nova21-14:43993 snodgras@nova.its.iastate.edu
```


*Split GVCF:*
```
for i in *.swap.maf; do
	n=$(echo $i | cut -f 2 -d "_") 
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/04.maf2gvcf.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.fasta ${n} >> ../scripts/maf2gvcf.cmds.txt
done
cd ../scripts/
#edit makeSLURM.py to be only 1 hour of wall time
python makeSLURM.py 1 maf2gvcf.cmds.txt
cd -
for i in ../scripts/maf2gvcf*.sub ; do sbatch $i ;done
```
using `bcftools`
```
ml bcftools
#bcftools view -t chr:from-to input.file > output.file

bcftools view -t chr1:54368399-59657020 Sb313_TdFL_anchorwave.swap.gvcf.gz> test_TdFL.gvcf.gz
#made a file, but didn't find any output

bcftools view -t chr17:14996000-15036000 Sb313_TdFL_anchorwave.swap.gvcf.gz> test2_TdFL.gvcf.gz
#made a file, but didn't find any output (these coordinates were taken from the MAF so it should be in the gvcf...)

bcftools view Sb313_TdFL_anchorwave.swap.gvcf.gz -t chr17:14996000-15036000 > test2_TdFL.gvcf.gz
#still didn't work

bcftools view Sb313_TdFL_anchorwave.swap.gvcf.gz -t chr17 > test2_TdFL.gvcf.gz
#still didn't work


```
#Swap coordinates of split MAF back to original





#Split the genomes based off parsed genespace results
reformat to bed format and split by consensusBin (last column)
```
for i in *_results.tsv ; do
	n=${i%_results.tsv}
	awk -v OFS='\t' '$8 ~ /bin1/ {print $1,$2,$3,$5}' $i > ${n}.bin1.bed
	awk -v OFS='\t' '$8 ~ /bin2/ {print $1,$2,$3,$5}' $i > ${n}.bin2.bed
done

for i in *.bed ; do 
	sort $i | uniq > uniq.${i} 
done
```
Then split the assembly files using the bed coordinates
```splitGenomes.sh
#!/bin/bash
outname=$1
assembly=$2
bed=$3

#module load bedtools2
bedtools getfasta -fo split_genome_assemblies/${outname} -fi ${assembly} -bed ${bed}
```

```
module load bedtools2 #this way it doesn't have to reload everytime
bedtools getfasta -fo split_genome_assemblies/TdFL.bin2.fasta -fi assemblies_final/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta -bed parsing_coordinates/uniq.zmB73vsTdFL.bin2.bed
bedtools getfasta -fo split_genome_assemblies/TdFL.bin1.fasta -fi assemblies_final/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta -bed parsing_coordinates/uniq.zmB73vsTdFL.bin1.bed 

bash scripts/splitGenomes.sh ZdGigi.bin1.fasta assemblies_final/Zd-Gigi-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZdGigi.bin1.bed
bash scripts/splitGenomes.sh ZdGigi.bin2.fasta assemblies_final/Zd-Gigi-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZdGigi.bin2.bed
bash scripts/splitGenomes.sh ZdMomo.bin1.fasta assemblies_final/Zd-Momo-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZdMomo.bin1.bed
bash scripts/splitGenomes.sh ZhRIMHU001.bin1.fasta assemblies_final/Zh-RIMHU001-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZhRIMHU001.bin1.bed
bash scripts/splitGenomes.sh ZnPI615697.bin1.fasta assemblies_final/Zn-PI615697-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZnPI615697.bin1.bed
bash scripts/splitGenomes.sh ZvTIL01.bin1.fasta assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZvTIL01.bin1.bed
bash scripts/splitGenomes.sh ZvTIL11.bin1.fasta assemblies_final/Zv-TIL11-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZvTIL11.bin1.bed
bash scripts/splitGenomes.sh ZxTIL18.bin1.fasta assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZxTIL18.bin1.bed
bash scripts/splitGenomes.sh ZxTIL25.bin1.fasta assemblies_final/Zx-TIL25-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZxTIL25.bin1.bed
bash scripts/splitGenomes.sh ZmB73.bin1.fasta NAM-assemblies/Zm-B73.fasta parsing_coordinates/uniq.zmB73vsZmB73.bin1.bed      
bash scripts/splitGenomes.sh ZmB97.bin1.fasta NAM-assemblies/Zm-B97.fasta parsing_coordinates/uniq.zmB73vsZmB97.bin1.bed      
bash scripts/splitGenomes.sh ZmCML103.bin1.fasta NAM-assemblies/Zm-CML103.fasta parsing_coordinates/uniq.zmB73vsZmCML103.bin1.bed   
bash scripts/splitGenomes.sh ZmCML228.bin1.fasta NAM-assemblies/Zm-CML228.fasta parsing_coordinates/uniq.zmB73vsZmCML228.bin1.bed    
bash scripts/splitGenomes.sh ZmCML247.bin1.fasta NAM-assemblies/Zm-CML247.fasta parsing_coordinates/uniq.zmB73vsZmCML247.bin1.bed    
bash scripts/splitGenomes.sh ZmCML277.bin1.fasta NAM-assemblies/Zm-CML277.fasta parsing_coordinates/uniq.zmB73vsZmCML277.bin1.bed    
bash scripts/splitGenomes.sh ZmCML322.bin1.fasta NAM-assemblies/Zm-CML322.fasta parsing_coordinates/uniq.zmB73vsZmCML322.bin1.bed   
bash scripts/splitGenomes.sh ZmCML333.bin1.fasta NAM-assemblies/Zm-CML333.fasta parsing_coordinates/uniq.zmB73vsZmCML333.bin1.bed    
bash scripts/splitGenomes.sh ZmCML52.bin1.fasta NAM-assemblies/Zm-CML52.fasta parsing_coordinates/uniq.zmB73vsZmCML52.bin1.bed 
bash scripts/splitGenomes.sh ZmCML69.bin1.fasta NAM-assemblies/Zm-CML69.fasta parsing_coordinates/uniq.zmB73vsZmCML69.bin1.bed
bash scripts/splitGenomes.sh ZmHP301.bin1.fasta NAM-assemblies/Zm-HP301.fasta parsing_coordinates/uniq.zmB73vsZmHP301.bin1.bed
bash scripts/splitGenomes.sh ZmIL14H.bin1.fasta NAM-assemblies/Zm-IL14H.fasta parsing_coordinates/uniq.zmB73vsZmIL14H.bin1.bed
bash scripts/splitGenomes.sh ZmKi11.bin1.fasta NAM-assemblies/Zm-Ki11.fasta parsing_coordinates/uniq.zmB73vsZmKi11.bin1.bed
bash scripts/splitGenomes.sh ZmKi3.bin1.fasta NAM-assemblies/Zm-Ki3.fasta parsing_coordinates/uniq.zmB73vsZmKi3.bin1.bed
bash scripts/splitGenomes.sh ZmKy21.bin1.fasta NAM-assemblies/Zm-Ky21.fasta parsing_coordinates/uniq.zmB73vsZmKy21.bin1.bed
bash scripts/splitGenomes.sh ZmM162W.bin1.fasta NAM-assemblies/Zm-M162W.fasta parsing_coordinates/uniq.zmB73vsZmM162W.bin1.bed
bash scripts/splitGenomes.sh ZmM37W.bin1.fasta NAM-assemblies/Zm-M37W.fasta parsing_coordinates/uniq.zmB73vsZmM37W.bin1.bed
bash scripts/splitGenomes.sh ZmMo18W.bin1.fasta NAM-assemblies/Zm-Mo18W.fasta parsing_coordinates/uniq.zmB73vsZmMo18W.bin1.bed
bash scripts/splitGenomes.sh ZmMS71.bin1.fasta NAM-assemblies/Zm-MS71.fasta parsing_coordinates/uniq.zmB73vsZmMS71.bin1.bed
bash scripts/splitGenomes.sh ZmNC350.bin1.fasta NAM-assemblies/Zm-NC350.fasta parsing_coordinates/uniq.zmB73vsZmNC350.bin1.bed 
bash scripts/splitGenomes.sh ZmNC358.bin1.fasta NAM-assemblies/Zm-NC358.fasta parsing_coordinates/uniq.zmB73vsZmNC358.bin1.bed
bash scripts/splitGenomes.sh ZmOh43.bin1.fasta NAM-assemblies/Zm-Oh43.fasta parsing_coordinates/uniq.zmB73vsZmOh43.bin1.bed
bash scripts/splitGenomes.sh ZmOh7b.bin1.fasta NAM-assemblies/Zm-Oh7b.fasta parsing_coordinates/uniq.zmB73vsZmOh7b.bin1.bed
bash scripts/splitGenomes.sh ZmP39.bin1.fasta NAM-assemblies/Zm-P39.fasta parsing_coordinates/uniq.zmB73vsZmP39.bin1.bed
bash scripts/splitGenomes.sh ZmTx303.bin1.fasta NAM-assemblies/Zm-Tx303.fasta parsing_coordinates/uniq.zmB73vsZmTx303.bin1.bed
bash scripts/splitGenomes.sh ZmTzi8.bin1.fasta NAM-assemblies/Zm-Tzi8.fasta parsing_coordinates/uniq.zmB73vsZmTzi8.bin1.bed

bash scripts/splitGenomes.sh ZdMomo.bin2.fasta assemblies_final/Zd-Momo-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZdMomo.bin2.bed
bash scripts/splitGenomes.sh ZhRIMHU001.bin2.fasta assemblies_final/Zh-RIMHU001-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZhRIMHU001.bin2.bed
bash scripts/splitGenomes.sh ZnPI615697.bin2.fasta assemblies_final/Zn-PI615697-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZnPI615697.bin2.bed
bash scripts/splitGenomes.sh ZvTIL01.bin2.fasta assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZvTIL01.bin2.bed
bash scripts/splitGenomes.sh ZvTIL11.bin2.fasta assemblies_final/Zv-TIL11-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZvTIL11.bin2.bed
bash scripts/splitGenomes.sh ZxTIL18.bin2.fasta assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZxTIL18.bin2.bed
bash scripts/splitGenomes.sh ZxTIL25.bin2.fasta assemblies_final/Zx-TIL25-REFERENCE-PanAnd-1.0.fasta parsing_coordinates/uniq.zmB73vsZxTIL25.bin2.bed
bash scripts/splitGenomes.sh ZmB73.bin2.fasta NAM-assemblies/Zm-B73.fasta parsing_coordinates/uniq.zmB73vsZmB73.bin2.bed      
bash scripts/splitGenomes.sh ZmB97.bin2.fasta NAM-assemblies/Zm-B97.fasta parsing_coordinates/uniq.zmB73vsZmB97.bin2.bed      
bash scripts/splitGenomes.sh ZmCML103.bin2.fasta NAM-assemblies/Zm-CML103.fasta parsing_coordinates/uniq.zmB73vsZmCML103.bin2.bed   
bash scripts/splitGenomes.sh ZmCML228.bin2.fasta NAM-assemblies/Zm-CML228.fasta parsing_coordinates/uniq.zmB73vsZmCML228.bin2.bed    
bash scripts/splitGenomes.sh ZmCML247.bin2.fasta NAM-assemblies/Zm-CML247.fasta parsing_coordinates/uniq.zmB73vsZmCML247.bin2.bed    
bash scripts/splitGenomes.sh ZmCML277.bin2.fasta NAM-assemblies/Zm-CML277.fasta parsing_coordinates/uniq.zmB73vsZmCML277.bin2.bed    
bash scripts/splitGenomes.sh ZmCML322.bin2.fasta NAM-assemblies/Zm-CML322.fasta parsing_coordinates/uniq.zmB73vsZmCML322.bin2.bed   
bash scripts/splitGenomes.sh ZmCML333.bin2.fasta NAM-assemblies/Zm-CML333.fasta parsing_coordinates/uniq.zmB73vsZmCML333.bin2.bed    
bash scripts/splitGenomes.sh ZmCML52.bin2.fasta NAM-assemblies/Zm-CML52.fasta parsing_coordinates/uniq.zmB73vsZmCML52.bin2.bed 
bash scripts/splitGenomes.sh ZmCML69.bin2.fasta NAM-assemblies/Zm-CML69.fasta parsing_coordinates/uniq.zmB73vsZmCML69.bin2.bed
bash scripts/splitGenomes.sh ZmHP301.bin2.fasta NAM-assemblies/Zm-HP301.fasta parsing_coordinates/uniq.zmB73vsZmHP301.bin2.bed
bash scripts/splitGenomes.sh ZmIL14H.bin2.fasta NAM-assemblies/Zm-IL14H.fasta parsing_coordinates/uniq.zmB73vsZmIL14H.bin2.bed
bash scripts/splitGenomes.sh ZmKi11.bin2.fasta NAM-assemblies/Zm-Ki11.fasta parsing_coordinates/uniq.zmB73vsZmKi11.bin2.bed
bash scripts/splitGenomes.sh ZmKi3.bin2.fasta NAM-assemblies/Zm-Ki3.fasta parsing_coordinates/uniq.zmB73vsZmKi3.bin2.bed
bash scripts/splitGenomes.sh ZmKy21.bin2.fasta NAM-assemblies/Zm-Ky21.fasta parsing_coordinates/uniq.zmB73vsZmKy21.bin2.bed
bash scripts/splitGenomes.sh ZmM162W.bin2.fasta NAM-assemblies/Zm-M162W.fasta parsing_coordinates/uniq.zmB73vsZmM162W.bin2.bed
bash scripts/splitGenomes.sh ZmM37W.bin2.fasta NAM-assemblies/Zm-M37W.fasta parsing_coordinates/uniq.zmB73vsZmM37W.bin2.bed
bash scripts/splitGenomes.sh ZmMo18W.bin2.fasta NAM-assemblies/Zm-Mo18W.fasta parsing_coordinates/uniq.zmB73vsZmMo18W.bin2.bed
bash scripts/splitGenomes.sh ZmMS71.bin2.fasta NAM-assemblies/Zm-MS71.fasta parsing_coordinates/uniq.zmB73vsZmMS71.bin2.bed
bash scripts/splitGenomes.sh ZmNC350.bin2.fasta NAM-assemblies/Zm-NC350.fasta parsing_coordinates/uniq.zmB73vsZmNC350.bin2.bed 
bash scripts/splitGenomes.sh ZmNC358.bin2.fasta NAM-assemblies/Zm-NC358.fasta parsing_coordinates/uniq.zmB73vsZmNC358.bin2.bed
bash scripts/splitGenomes.sh ZmOh43.bin2.fasta NAM-assemblies/Zm-Oh43.fasta parsing_coordinates/uniq.zmB73vsZmOh43.bin2.bed
bash scripts/splitGenomes.sh ZmOh7b.bin2.fasta NAM-assemblies/Zm-Oh7b.fasta parsing_coordinates/uniq.zmB73vsZmOh7b.bin2.bed
bash scripts/splitGenomes.sh ZmP39.bin2.fasta NAM-assemblies/Zm-P39.fasta parsing_coordinates/uniq.zmB73vsZmP39.bin2.bed
bash scripts/splitGenomes.sh ZmTx303.bin2.fasta NAM-assemblies/Zm-Tx303.fasta parsing_coordinates/uniq.zmB73vsZmTx303.bin2.bed
bash scripts/splitGenomes.sh ZmTzi8.bin2.fasta NAM-assemblies/Zm-Tzi8.fasta parsing_coordinates/uniq.zmB73vsZmTzi8.bin2.bed

```
Double check that there weren't any typos above

```
#check the md5sum and make sure there are as many unique sums as there are files
md5sum *.fasta | cut -f 1 -d " " | sort | uniq | wc -l

ls *.fasta | wc -l 
```
They should both be 70!

#Running AnchorWave on the split query genomes
AnchorWave requires:
- Reference (Sorghum) assembly fasta
- Reference (Sorghum) gene gff3
- Query (Tripsacinae) assembly fasta 

##Step 1: making the reference CDS files
Script: `slurm_01.makeRefCDS.sh`
Estimated Time to Run: _1 minute_

```
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave gff2seq \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \ #ref genome sequence
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \ #ref genome annotation as gff or gff3
	-o Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa #output file of longest CDS/exon for each gene
```

##Step 2: mapping with minimap2 to get anchorpoints
Script: `slurm_02.minimapAnchorPoints.sh` and `02.minimapAnchorPoints.sh` (the if statement has weird colors in vi, not sure if it'll run...)
Time to Run: 00:01:29:00 ish

Make a directory for the sam files to go to: 
```
mkdir sam_files
```
And make a softlink for the Av genome to the split genomes file
```
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta split_genome_assemblies/Av.fasta
```
```bash script
FILE1=$1 #Genomic fasta file ID for maize genome (genomeID)) (example: Zm-B73-REFERENCE-NAM-5.0); be sure to make certain the fasta file extension in the script matches user file extension i.e. '.fasta' vs '.fa', etc)
FileName=${FILE1#split_genome_assemblies/}
ml minimap2

#Do I need to redo the minimap to sorghum to itself everytime? Could cut down with a if file exists exemption
if [[ ! -f "sam_files/Sb313.ref.sam" ]] ; then
	minimap2 \
		-x splice \
		-t 11 \
		-k 12 \
		-a \
		-p 0.4 \
		-N 20 \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > sam_files/Sb313.ref.sam
fi

minimap2 \
	-x splice \ #Applies multiple options at once: long-read spliced alignment; long deletions are taken as introns and represented as "n" in CIGAR, long insertions are disabled, indel gap costs are different during chaining; computing 'ms' tag ignores introns to demote hits to pseudogenes
	-t 11 \ #threads
	-k 12 \ #minimizer k-mer length
	-a \ #generate CIGAR and ouput alignmnets in SAM format
	-p 0.4 \ #Minimal secondary-to-primary score ratio to output secondary mappings
	-N 20 \ #Output at most INT secondary alignments
	${FILE1} \ #target fasta
	Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > sam_files/${FileName%.fasta}_Sb313.cds.sam #query fasta and output file
```

```slurm script
for i in split_genome_assemblies/*.fasta ; do
        bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/02.minimapAnchorPoints.sh ${i}
;done
```

##Step 3: Anchorwave Alignment
Running Anchorwave will take ~8 hours without -f and ~24 hours with -f
Needs ~300G of memory for 5 threads
Will not run on gzipped fasta files
Timed out after 12 hours
-r = ref genome (-R is max ref alignment coverage)
-s = target genome sequence (-Q = query max alignment coverage 


```scripts/03.AnchorWaveAlignment.sh
FILE1=$1
FileName=${FILE1#split_genome_assemblies/}
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \ 
	-as Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \ 
	-a sam_files/${FileName}_Sb313.cds.sam \ 
	-ar sam_files/Sb313.ref.sam \ 
	-s ${FILE1}.fasta \ 
	-n AnchorWave_output/Sb313_${FileName}_anchorwave.anchors \ 
	-R 1 \ 
	-Q 1 \
	-o AnchorWave_output/Sb313_${FileName}_anchorwave.maf \ 
	-t 5 > AnchorWave_logs/Sb313_${FileName}_anchorwave.log 2>&1 
####
#slurm job
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="AnchorWaveAlignment"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive

bash scripts/03.AnchorWaveAlignment.sh assemblies_final/Zv-TIL01-Reference-PanAnd-2.0
```
To make multiple slurm scripts for parallel job submission
```
for i in split_genome_assemblies/*.fasta ; do
	echo bash scripts/03.AnchorWaveAlignment.sh ${i%.fasta} >> scripts/AnchorWave.cmds.txt
done
python makeSLURM.py 1 AnchorWave.cmds.txt
```

```
#but to set the other ones up to run if it completes ok
for i in {1..70}; do
	sbatch --dependency=afterok:4695887 scripts/AnchorWave.cmds_${i}.sub ;
done
```
***RAN OUT OF ROOM IN WORK/LAS/ SO I'M RUNNING EVERY JOB AFTER 22 IN PTMP*
```
cd /ptmp/LAS/snodgras/
ln -s /work/LAS/mhufford-lab/snodgras/Fractionation/ .

for i in {23..70} ; do sbatch scripts/AnchorWave.cmds_${i}.sub ; done
```

##Step 4: convert MAF to GVCF formats
Will need anchorwave and tassel to run
scripts below: `04.maf2gvcf.sh`
It also requires the phg singularity
```
#phg singularity installation
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load singularity
singularity pull --name phg_latest.sif docker://maizegenetics/phg
```
```
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
singularity exec phg_latest.sif /tassel-5-standalone/run_pipeline.pl \
	-Xmx300g \
	-MAFToGVCFPlugin \
	-referenceFasta ${REFfasta} \
	-mafFile Sb313_${REFname}_anchorwave.maf \
	-sampleName Sb313_${REFname} \
	-gvcfOutput Sb313_${REFname}_anchorwave.gvcf \
	-fillGaps true

##The command "-Xmx300g" demands that Java has enough RAM for the job to be run; in this case, 300g

###RUN IN /AnchorWave_output/

#SLURM COMMAND
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="MAF2GVCF"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash /ptmp/LAS/snodgras/Fractionation/scripts/04.maf2gvcf.sh /ptmp/LAS/snodgras/Fractionation/Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa Av 
```
`sbatch ../scripts/slurm_04.maf2gvcf.sh`

To make it run parallel
```
for i in *bin*.maf ; do
	n=$(echo $i | cut -f 2 -d "_") 
	echo bash /ptmp/LAS/snodgras/Fractionation/scripts/04.maf2gvcf.sh /ptmp/LAS/snodgras/Fractionation/Sb313.fasta ${n} >> ../scripts/maf2gvcf.cmds.txt
done
cd ../scripts/
#edit makeSLURM.py to be only 1 hour of wall time
python makeSLURM.py 1 maf2gvcf.cmds.txt
cd -
for i in ../scripts/maf2gvcf*.sub ; do sbatch $i ;done
```

I had to copy the sb313 fasta again because it was saying it couldn't find it when running the slurm

##Step 5: Trimming GVCF with GATK
Original code from Elena Jiang, altered by Samantha Snodgrass
```slurm_05.setupRef.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=1:00:00
#SBATCH --job-name=gatk-trim
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

#load module
module load samtools
module load gatk

sed -e 's/Chr//g' /ptmp/LAS/snodgras/Fractionation/Sb313.fasta > Sb313.clean.fasta
samtools faidx Sb313.clean.fasta
gatk CreateSequenceDictionary -R Sb313.clean.fasta
```
Change the `makeSLURM.py` to match the above slurm settings
```
for i in *gvcf.gz ; do 
echo module load samtools \; module load gatk \; gatk LeftAlignAndTrimVariants -O trimmed_${i} -V ${i} -R /ptmp/LAS/snodgras/Fractionation/Sb313.clean.fasta --max-indel-length 9101263 &> gf_stdout.txt >> ../scripts/trimGVCF.cmds.txt
done
cd ../scripts/
#edit makeSLURM.py to be only 1 hour of wall time
python makeSLURM.py 1 trimGVCF.cmds.txt
cd -
sbatch ../scripts/trimGVCF.cmds_0.sub
for i in {1..70} ; do sbatch --dependency=afterok: ../scripts/trimGVCF.cmds_${i}.sub ;done
```

##Step 6: gatk gvcf to vcf
_Code from Elena Jiang_
You can also reference baoxing’s script here https://github.com/baoxingsong/AnchorWave/blob/master/doc/GATK.md#14-merge-gvcf-files-into-a-gatk-database

*Split the reference genome into separate chromosome fastas*
In this case Sb313
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < Sb313.clean.fasta > Sb313.linear.clean.fasta
#manually remove empty first line
grep -A1 "^>01" Sb313.linear.clean.fasta > Sb313.chr1.fasta
grep -A1 "^>02" Sb313.linear.clean.fasta > Sb313.chr2.fasta
grep -A1 "^>03" Sb313.linear.clean.fasta > Sb313.chr3.fasta
grep -A1 "^>04" Sb313.linear.clean.fasta > Sb313.chr4.fasta
grep -A1 "^>05" Sb313.linear.clean.fasta > Sb313.chr5.fasta
grep -A1 "^>06" Sb313.linear.clean.fasta > Sb313.chr6.fasta
grep -A1 "^>07" Sb313.linear.clean.fasta > Sb313.chr7.fasta
grep -A1 "^>08" Sb313.linear.clean.fasta > Sb313.chr8.fasta
grep -A1 "^>09" Sb313.linear.clean.fasta > Sb313.chr9.fasta
grep -A1 "^>10" Sb313.linear.clean.fasta > Sb313.chr10.fasta
mv Sb313.chr*.fasta AnchorWave_Output/.
```

*-V for every trimmed gvcf file created from the previous step (#5)*
*-L is the number of the chromosome (chromosome name) that is being used*
```slurm_06.01.gvcf2vcf.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=96:00:00
#SBATCH --job-name=gatk-chr1
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml gatk
ref=Sb313.chr1.fasta
ml samtools
# index
samtools faidx $ref
# dict
gatk CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${ref%.*}.dict
# db import
gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
-V trimmed_Sb313_Av_anchorwave.gvcf.gz \
-V trimmed_Sb313_TdFL.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_TdFL.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZdGigi.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZdGigi.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZdMomo.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZdMomo.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZhRIMHU001.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZhRIMHU001.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmB73.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmB73.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmB97.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmB97.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML103.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML103.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML228.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML228.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML247.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML247.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML277.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML277.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML322.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML322.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML333.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML333.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML52.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML52.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML69.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmCML69.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmHP301.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmHP301.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmIL14H.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmIL14H.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmKi11.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmKi11.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmKi3.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmKi3.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmKy21.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmKy21.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmM162W.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmM162W.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmM37W.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmM37W.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmMo18W.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmMo18W.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmMS71.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmMS71.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmNC350.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmNC350.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmNC358.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmNC358.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmOh43.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmOh43.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmOh7b.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmOh7b.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmP39.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmP39.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmTx303.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmTx303.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmTzi8.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZmTzi8.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZnPI615697.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZnPI615697.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZvTIL01.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZvTIL01.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZvTIL11.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZvTIL11.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZxTIL18.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZxTIL18.bin2_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZxTIL25.bin1_anchorwave.gvcf.gz \
-V trimmed_Sb313_ZxTIL25.bin2_anchorwave.gvcf.gz \
--batch-size 1 \
--genomicsdb-workspace-path chr1_gatkDBimport \
-L 1 \
--genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 
#--tmp-dir $TMPDIR

# vcf output
gatk --java-options "-Xmx50g" GenotypeGVCFs \
-R $ref \
-V gendb://chr1_gatkDBimport \
--cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000 \
-O chr1.vcf
```
*ERRORS*
For Chr01, it said some of the -V arguments were not positional
- checked for extra spaces after the `\`
- checked to make sure that there were Chr01 coordinates in the maf
```
for i in ZdGigi.bin2 ZmCML103.bin2 ZmCML333.bin1 ZmKi3.bin1 ; do grep "Chr01" Sb313_${i}_anchorwave.maf | wc -l ; done
0
0
11
6
```
- Try the script with only 3 -V (Av and the two bins for TdFL)

Make one slurm script for each chromosome, changing the 
- reference fasta name
- L option
- workspace path
- vcf V option
- vcf O option

```
for i in 02 03 04 05 06 07 08 09 10 ; do cp slurm_06.01.gvcf2vcf.sh slurm_06.${i}.gvcf2vcf.sh; done
sed -i 's/chr1/chr2/g' slurm_06.02.gvcf2vcf.sh
sed -i 's/chr1/chr3/g' slurm_06.03.gvcf2vcf.sh
sed -i 's/chr1/chr4/g' slurm_06.04.gvcf2vcf.sh
sed -i 's/chr1/chr5/g' slurm_06.05.gvcf2vcf.sh
sed -i 's/chr1/chr6/g' slurm_06.06.gvcf2vcf.sh
sed -i 's/chr1/chr7/g' slurm_06.07.gvcf2vcf.sh
sed -i 's/chr1/chr8/g' slurm_06.08.gvcf2vcf.sh
sed -i 's/chr1/chr9/g' slurm_06.09.gvcf2vcf.sh
sed -i 's/chr1/chr10/g' slurm_06.10.gvcf2vcf.sh

#manually change -L
```



#####OLD CODE#####

# ANCHORWAVE PIPELINE (used)
###Step 0: getting the files into the right directory

*Sorghum Genome Files*
This is where the sorghum genome and annotations come from: 
`https://data.jgi.doe.gov/refine-download/phytozome?organism=Sbicolor&expanded=Phytozome-313`
I'll need
```
Sbicolor_313_v3.0.fa.gz
Sbicolor_313_v3.1.gene.gff3
```
From Maggie: 
```
fasta file: Sbicolor_313_v3.0
gff file: Sbicolor_313_v3.1

from Phytozome V12
```

*Installing AnchorWave container*
```
ssh novadtn
cd /work/LAS/mhufford-lab/snodgras/Fractionation/
wget -O anchorwave.sif https://iastate.box.com/shared/static/2wy98iuvafpb39s98bjg1ljrr21nejk6.sif
```
Note: whenever I use AnchorWave, I'll have to do it like this:
```
ml singularity
singularity exec --bind $PWD anchorwave.sif your_command options -arguments etc
```
`$PWD` is the path with the input files; without it (or the wrong path), the singularity won't be able to "see" our input files
See notes in `Anchorwave_protocol_with_sorghum_as_reference.txt`

##Step 1: making the reference CDS files
Script: `slurm_01.makeRefCDS.sh`
Estimated Time to Run: _1 minute_
```
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave gff2seq \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \ #ref genome sequence
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \ #ref genome annotation as gff or gff3
	-o Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa #output file of longest CDS/exon for each gene
```

##Step 2: mapping with minimap2 to get anchorpoints
Script: `slurm_02.minimapAnchorPoints.sh` and `02.minimapAnchorPoints.sh` (the if statement has weird colors in vi, not sure if it'll run...)
Estimated Time to Run: _~35 minutes if unzipping all files_
For running on Panand use directory `assemblies_final/`
For running on NAM use directory `NAM-assemblies/` _took ~50 minutes with  unzipped files_
```
FILE1=$1 #Genomic fasta file ID for maize genome (genomeID)) (example: Zm-B73-REFERENCE-NAM-5.0); be sure to make certain the fasta file extension in the script matches user file extension i.e. '.fasta' vs '.fa', etc)
FileName=${FILE1#assemblies_final/}
ml minimap2

minimap2 \
	-x splice \ #Applies multiple options at once: long-read spliced alignment; long deletions are taken as introns and represented as "n" in CIGAR, long insertions are disabled, indel gap costs are different during chaining; computing 'ms' tag ignores introns to demote hits to pseudogenes
	-t 11 \ #threads
	-k 12 \ #minimizer k-mer length
	-a \ #generate CIGAR and ouput alignmnets in SAM format
	-p 0.4 \ #Minimal secondary-to-primary score ratio to output secondary mappings
	-N 20 \ #Output at most INT secondary alignments
	${FILE1}.fasta \ #target fasta
	Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > ${FileName}_Sbicolor_313_v3.0_cds.sam #query fasta and output file
#Do I need to redo the minimap to sorghum to itself everytime? Could cut down with a if file exists exemption
if [[ ! -f "Sbicolor_313_v3.0_ref.sam" ]] ; then
	minimap2 \
		-x splice \
		-t 11 \
		-k 12 \
		-a \
		-p 0.4 \
		-N 20 \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
		Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa > Sbicolor_313_v3.0_ref.sam
fi
```
Trying it on `panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta.gz`
It worked, so doing it in a loop for the other genomes:
```
for i in assemblies_final/*.gz ; do
        gunzip $i
        bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/02.minimapAnchorPoints.sh ${i%.fasta.gz}
;done
```
Got this warning
```
[WARNING] SAM output is malformated due to internal @SQ lines. Please add option --no-sam-sq or filter afterwards.
```
Plot full length CDS mapping results? 
Need the `alignmnetToDotplot.pl` script from anchorwave github
```
#perl AlignmentToDotplot.pl Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 query.sam > query.tab
perl scripts/AlignmentToDotplot.pl Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 Zv-TIL01-Reference-PanAnd-2.0_Sbicolor_313_v3.0_cds.sam > Zv-TIL01.Sbicolor.tab
```
Get chromosome info: 
```
samtools faidx Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa
samtools faidx panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta.bgz #had to gunzip -c panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta.gz | bgzip > panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta.bgz
```
Then use R to plot
```
library(ggplot2)
library(compiler)
enableJIT(3)
library(ggplot2)
library("Cairo")
changetoM <- function ( position ){
	position=position/1000000; 
	paste(position, "M", sep="")
}
data=read.table("query.tab")
data = data[which(data$V1 %in% c(REF CHROMOSOME IDS AS CHARACTER STRINGS)),]
data = data[which(data$V3 %in% c(QUERY CHROMOSOME IDS AS CHARACTER STRINGS)),]
data$V1 = factor(data$V1, levels=c(REF CHROMOSOME IDS IN ORDER))
data$V3 = factor(data$V3, levels=c(QUERY CHROMOSOME IDS IN ORDER))
ggplot(data = data, aes(x=V4, y=V2))+geom_point(size=0.5, aes(color=V5))+facet_grid(V1~V3, scales="free", space="free") + theme_grey(base_size=120)+
	labs(x="QUERY NAME", y="REF NAME")+scale_x_continuous(labels=changetoM)+scale_y_continuous(labels=changetoM)+
	theme(axis.line = element_blank(),
		panel.background = element_blank(),
		panel.border = element_rect(fill=NA,color="black",size=0.5,linetype="solid"),
		axis.text.y = element_text(colour = "black"),
		legend.position="none",
		axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
#Somehow return the plot as an output file
```
This was a weird dotplot so trying to using 1.2.5 in the walkthrough that uses the anchors
This should take ~15 minutes if all the other files have been made `slurm_02.5.AnchorDotplot.sh`
(Only creating the anchor files; if anchor files are already created, then use those)
```
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali -i Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \
	-r Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \
	-a Zv-TIL01-Reference-PanAnd-2.0_Sbicolor_313_v3.0_cds.sam \
	-as Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-ar Sbicolor_313_v3.0_ref.sam \
	-s panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta.gz \
	-n align1.anchors \
	-R 2 -Q 1 -ns
```
Then use R like above

##Step 3: Anchorwave Alignment
Running Anchorwave will take ~8 hours without -f and ~24 hours with -f
Needs ~300G of memory for 5 threads
Will not run on gzipped fasta files
Timed out after 12 hours
-r = ref genome (-R is max ref alignment coverage)
-s = target genome sequence (-Q = query max alignment coverage 

For running on Panand use directory `assemblies_final/`
For running on NAM use directory `NAM-assemblies/`

```
FILE1=$1
FileName=${FILE1#assemblies_final/}
ml singularity
singularity exec --bind $PWD anchorwave.sif anchorwave proali \
	-i Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 \ 
	-as Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0_cds.fa \
	-r Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.0.fa \ 
	-a ${FileName}_Sbicolor_313_v3.0_cds.sam \ 
	-ar Sbicolor_313_v3.0_ref.sam \ 
	-s ${FILE1}.fasta \ 
	-n Sbicolor_313_v3.0_${FileName}_anchorwave.anchors \ 
	-R 2 \ 
	-Q 1 \
	-o Sbicolor_313_v3.0_${FileName}_anchorwave.maf \ 
	-t 5 > Sbicolor_313_v3.0_${FileName}_anchorwave.log 2>&1 
####
#slurm job
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="AnchorWaveAlignment"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --exclusive

bash scripts/03.AnchorWaveAlignment.sh assemblies_final/Zv-TIL01-Reference-PanAnd-2.0
```
To make multiple slurm scripts for parallel job submission
```
for i in assemblies_final/*.fasta ; do
	echo bash scripts/03.AnchorWaveAlignment.sh ${i%.fasta} >> scripts/AnchorWave.cmds.txt
done
python makeSLURM.py 1 AnchorWave.cmds.txt
```
Testing with sub command 9 (Zx TIL25)
```
#but to set the other ones up to run if it completes ok
for i in {0..8}; do
	sbatch --dependency=afterok:4196500 scripts/AnchorWave.cmds_${i}.sub ;
done
```

##Step 4: convert MAF to GVCF formats
Will need anchorwave and tassel to run
scripts from AnchorWave: `anchorwave-maf-swap.py`
scripts below: `04.maf2gvcf.sh`
It also requires the phg singularity
```
#phg singularity installation
ml singularity
export SINGULARITY_CACHEDIR=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR
singularity pull docker://maizegenetics/phg
```
```
REFfasta=$1 #Maize genome fasta
REFname=$2 #name without extensions

##Anchorwave swap (this takes only a few minutes):

cat Sbicolor_313_v3.0_${REFname}_anchorwave.maf | python2 scripts/anchorwave-maf-swap.py > Sbicolor_313_v3.0_${REFname}_anchorwave.swap.maf

##Convert maf to gvcf (this takes about 15 or so minutes):

#MAFToGVCFPlugin <options>
#-referenceFasta <Reference Fasta> : Input Reference Fasta (required)
#-mafFile <Maf File> : Input MAF file.  Please note that this needs to be a MAF file with 2 samples.  The first will be assumed to be the Reference and the second will be the assembly. (required)
#-sampleName <Sample Name> : Sample Name to write to the GVCF file as the genome header or ID (required)
#-gvcfOutput <Gvcf Output> : Output GVCF file name (required)
#-fillGaps <true | false> : When true, if the maf file does not fully cover the reference genome any gaps in coverage will be #filled in with reference blocks. This is necessary if the resulting GVCFs are to be combined. (Default: false)

module load singularity
singularity exec phg_latest.sif /tassel-5-standalone/run_pipeline.pl \
	-Xmx300g \
	-MAFToGVCFPlugin \
	-referenceFasta ${REFfasta} \
	-mafFile Sbicolor_313_v3.0_${REFname}_anchorwave.swap.maf \
	-sampleName Sbicolor_313_v3.0_${REFname} \
	-gvcfOutput Sbicolor_313_v3.0_${REFname}_anchorwave.swap.gvcf \
	-fillGaps true

##The command "-Xmx300g" demands that Java has enough RAM for the job to be run; in this case, 300g

#SLURM COMMAND
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=369G   # maximum memory per node
#SBATCH --job-name="MAF2GVCF"
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/04.maf2gvcf.sh panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta Zv-TIL01-Reference-PanAnd-2.0
```
To make it run parallel
```
for i in Zl-RIL003-Reference-PanAnd-2.0 Zv-TIL11-Reference-PanAnd-2.0 Zx-TIL18-Reference-PanAnd-2.0 Zx-TIL25-Reference-PanAnd-2.0 ; do
	echo bash scripts/04.maf2gvcf.sh panand-assemblies/${i}.fasta ${i} >> maf2gvcf.cmds.txt
done
python makeSLURM.py 1 maf2gvcf.cmds.txt
for i in scripts/maf2gvcf*.sub ; do sbatch $i ;done
```

# Parsing gvcf files
Currently, only have the `.swap.gvcf` so ref is the tripsacinae genome and query is the sorghum.
*Subsetting so only have deletions*
script: `05.a.cleaningGVCF.sh`; run `for i in *swap.gvcf ; do bash scripts/05.a.cleaningGVCF.sh $i ;done`
Takes ~13 minutes to run
```
gvcf=$1

grep "^#" $gvcf > clean.$gvcf
grep -v "^#" $gvcf | awk -v OFS='\t' ' $5 ~ /,/ {print $0}' - >> clean.$gvcf
sed -i 's/,<NON_REF>//g' clean.$gvcf
```
script: `05.b.delOnlyGVCF.sh`; run `for i in clean*.gvcf ; do bash scripts/05.b.delOnlyGVCF.sh $i ; done`
Takes ~3 minutes to run
```
gvcf=$1

grep "^#CHROM" $gvcf > delOnly.$gvcf
awk -v OFS='\t' '{if((length($4) == 1) && (length($5) >= 2)) print $0}' $gvcf >> delOnly.$gvcf 
```
`mv delOnly* swap_gvcf_parsing/.`
*trying to make unique identifiers to find shared deletions*
script: `05.c.splitSortGVCF.sh` `for i in swap_gvcf_parsing/delOnly*.gvcf; do bash scripts/05.c.splitSortGVCF.sh $i ; done`
Takes ~10 minutes to run
```
gvcf=$1

head -n 1 $gvcf | awk -v OFS='\t' '{print $0,"ASM_Chr","ASM_End","ASM_Start","ASM_Strand","END"}' - > ${gvcf%-Reference-PanAnd-2.0_anchorwave.swap.gvcf}.id.tmp

tail -n +2 $gvcf | 
	awk -v OFS='\t' '{split($8,a,/;/); split(a[1],b,/=/); split(a[2],c,/=/); split(a[3],d,/=/); split(a[4],e,/=/); split(a[5],f,/=/); print $0, b[2],c[2],d[2],e[2],f[2]}' - |
	sort -k11,13 >> ${gvcf%-Reference-PanAnd-2.0_anchorwave.swap.gvcf}.id.tmp
```
Check to see if there's any where multiple deletions are overlapping the same sorghum (ASM) coordinates
We do expect to have 2 unique tripsacinae pieces matching to each sorghum coordinate...
...but since we're looking at the deletions and not the present sequence...
...We still may have duplicate sorghum coordinates (if both copies were lost) or if deletion calls were overlapping (though shouldn't because of the way gvcf work...)
```
wc -l *.id.tmp
for i in *.id.tmp ; do cut -f 11-16 $i | uniq | wc -l ; done 
```
Found 
```
2559689 delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.id.tmp
   2554704 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.id.tmp
   2528519 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.id.tmp
   2530139 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.id.tmp
   2530063 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.id.tmp
  12703114 total

2552441 delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.id.tmp (7248 less rows)
2546750 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.id.tmp (7954 less rows)
2519477 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.id.tmp (9042 less rows)
2519364 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.id.tmp (10775 less rows)
2522263 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.id.tmp (7800 less rows)
12660295 total (42819 less rows)
```
let's try no cut (see if it's the same Sb coords, but diff Tripsacinae coords)
```
for i in *.id.tmp ; do echo $i; uniq $i | wc -l ; done

2559689 delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.id.tmp
2554704 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.id.tmp
2528519 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.id.tmp
2530139 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.id.tmp
2530063 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.id.tmp
```
Yep! This suggests to me that there are some either complete or partial losses shared across Sorghum coordinates.
That would be why the same Sb coords for a deletion event would map multiple times with unique Tripsacinae coords
Let's make sure that it's not something in the other fields causing the uniqueness:
```
for i in *.id.tmp ; do cut -f 1-2,11-15 $i | uniq |wc -l ; echo $i ;done

2559689 delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.id.tmp
2554704 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.id.tmp
2528519 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.id.tmp
2530139 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.id.tmp
2530063 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.id.tmp
```

script:`05.d.addIDs.sh` run as `for i in delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.id.tmp delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.id.tmp delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.id.tmp delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.id.tmp delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.id.tmp; do bash scripts/05.d.addIDs.sh swap_gvcf_parsing $i ; done`
Takes ~1 minute time to run
```
directory=$1
gvcf=$2

cd $directory

awk -v OFS='\t' '{print $11"."$13"."$12}' $gvcf > ${gvcf%.id.tmp}.id.list
awk -v OFS='\t' '{print $0,$11"."$13"."$12}' $gvcf > ${gvcf%.id.tmp}.id.full
```
```
cat *.id.list | sort | uniq  > Sb.coord.ID.key
cat *.id.list | sort | uniq | wc -l
9102996
wc -l *.id.list
  2559689 delOnly.clean.Sbicolor_313_v3.0_Zl-RIL003.id.list
  2554704 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL01.id.list
  2528519 delOnly.clean.Sbicolor_313_v3.0_Zv-TIL11.id.list
  2530139 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL18.id.list
  2530063 delOnly.clean.Sbicolor_313_v3.0_Zx-TIL25.id.list
 12703114 total
```
So once the duplicate deletions in the individual genomes are taken into account (12,660,295 total vs. 12,703,114)...
...there are 9,102,996 deletions with unique Sb coords across these 5 genomes...
...so there are 3,557,299 deletions with shared Sb coords across these 5 genomes.

To test out my R script to make the matrix:
```
for i in *.full ; do head -n 10000 $i > ${i%.id.full}.test ; done
```

## Redo pipeline to be more streamlined with less intermediate files
step 1: filter for only deletions and remove `<NON_REF>` character strings
step 2: split out the assembly coordinates and genotype information (Sorghum coordinates) and reformat to bed file with sorghum as the first columns
step 3: filter for deletions within sorghum exons
step 4: add IDs

Script: `05.A.cleanDelsOnly.sh`: `for i in *swap.gvcf; do bash scripts/05.A.cleanDelsOnly.sh $i ; done`
Adds header
Grabs everything non-header | 
prints only if there's a comma in ALT | 
removes `<NON_REF>` | 
prints only if the REF is length 1 and ALT is length 2+ > writes to "cleanDELs.GVCFFILE"
```
gvcf=$1

grep "^#" $gvcf | tail -n 1 > cleanDELs.$gvcf
grep -v "^#" $gvcf | awk -v OFS='\t' ' $5 ~ /,/ {print $0}' - | \
	sed 's/,<NON_REF>//g' - | \
	awk -v OFS='\t' '{if((length($4) == 1) && (length($5) >= 2)) print $0}' - >> cleanDELs.$gvcf
```

Script: `05.B.GVCF2BED.sh`: `for i in cleanDELs* ; do bash scripts/05.B.GVCF2BED.sh $i ; done`
Adds header
Takes everything that's not header and splits the 8th column with all the information and then sorts it
```
gvcf=$1

head -n 1 $gvcf | awk -v OFS='\t' '{print "ASM_Chr","ASM_Start","ASM_End","ASM_Strand",$0}' - > ${gvcf%-Reference-PanAnd-2.0_anchorwave.swap.gvcf}.bed

tail -n +2 $gvcf | 
	awk -v OFS='\t' '{split($8,a,/;/); split(a[1],b,/=/); split(a[2],c,/=/); split(a[3],d,/=/); split(a[4],e,/=/); split(a[5],f,/=/); print b[2],d[2],c[2],e[2],f[2],$0}' - |
	sort -k1,2 >> ${gvcf%-Reference-PanAnd-2.0_anchorwave.swap.gvcf}.bed
```

Script: `05.C.GenicAnchorsOnly.sh`
```
bed=$1
exon=$2
anchor=$3

#gets rid of all the lines starting with # and the header, then deals with any inversions in the ref coordinates
grep -v "^#" $anchor | tail -n +2 | awk -v OFS='\t' '{if ($2> $3) print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10; else print $0}' - > ${anchor%-Reference-PanAnd-2.0_anchorwave.anchors}.anchors.bed

ml bedtools2

#Add "Chr" to ASM_Chr cleanDEL.bed and rm columns ID, QUAL, FILTER, INFO
awk -v OFS='\t' '{print "Chr"$1,$2,$3,$4,$5,$6,$8,$9,$13,$14}' $bed | \
	tail -n +2 | \
	bedtools intersect -a - -b $exon -wa -wb | \
	bedtools intersect -a - -b ${anchor%-Reference-PanAnd-2.0_anchorwave.anchors}.anchors.bed -wa -wb | \
	cut -f 1-10,14,23 > gene.anchored.${bed}

#columns should be ASM_Chr, ASM_Start, ASM_End, ASM_Strand, #CHROM, POS, REF, ALT, FORMAT, TRIP-NAME, 
#exonCHR, exonStart, exonStop, SorgGeneID, Strand
#refChr	referenceStart	referenceEnd	queryChr	queryStart	queryEnd	strand	gene	blockIndex	score

#cut to keep the alignment coordinates and information (1-10), the gene/exon id (14), and the anchor type of the anchor (23)

```

Format of longest transcript bed file:
Chr01	start	stop	ID	strand

Format of anchor file: 
refChr	referenceStart	referenceEnd	queryChr	queryStart	queryEnd	strand	gene	blockIndex	score
Chr02	67783780	67786305	chr7	193243607	193246024	+	Sobic.002G301000.1.v3.1	1	0.982667
*Will need to get rid of lines starting with #*

Format of cleanDEL bed file:
ASM_Chr	ASM_Start	ASM_End	ASM_Strand	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO FORMAT	Zl-RIL003-Reference-PanAnd-2.0
01	10000077	10000078	-		1	371899900	.	A	AT	.	.	ASM_Chr=01;ASM_End=10000078;ASM_Start=10000077;ASM_Strand=-	GT:AD:DP	1:0,1,0:1
*ASM is actually the sorghum coordinates and CHROM POS is the Tripsacinae coordinates*

#Finding orthology groups to split MAF files before the GVCF conversion
To run jupyter notebooks on nova: https://nova-ood.its.iastate.edu/pun/sys/dashboard/
Help pages (https://researchit.las.iastate.edu/how-use-jupyterlab-open-ondemand)