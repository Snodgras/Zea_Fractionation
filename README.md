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

# ANCHORWAVE PIPELINE (used)
###Step 0: getting the files into the right directory
Have to do this on the dtn node to access LSS
```
rsync -rts --progress /lss/research/mhufford-lab/arnstrm/fractionation /ptmp/LAS/snodgras
```
Inspecting the files:
```
gunzip B97.gvcf.gz
```
*Teosinte genomes from PanAnd**
Zl = Zea luxurians
Zv = Zea mays spp. parviglumis
Zx = Zea mays spp. mexicana
Version meanings: <1 = contigs, 1 - < 2 = scaffolds, >= 2 = pseudomolecules
Highest number is most updated assembly

Dactyloides, luxurians, mexicana, and parv == pseudomolecules
Diploperennis == mess
Nicaraguensis, huehuetengensis, zopolitense == contig (would need to scaffold them)

These are the genome files names

*Copy assemblies to my local directory*
```
ssh novadtn
rsync -rts --progress /ptmp/arnstrm/LATEST_ANDRO_ASSEMBLIES/Z*2.0.fasta.gz /work/LAS/mhufford-lab/snodgras/Fractionation/panand-assemblies/.
```
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

##Running pSONIC
I've put each anchorwave set of files into specific directories for each of the genomes
`https://github.com/conJUSTover/pSONIC`

```
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanX
make
git clone https://github.com/conJUSTover/pSONIC
scp /work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2/stats/gff3/gff/Tripsacum* /work/LAS/mhufford-lab/snodgras/Fractionation/panand-annotations/.
scp /work/LAS/mhufford-lab/arnstrm/PanAnd/current-versions-genomes.v2/stats/gff3/gff/Zea* /work/LAS/mhufford-lab/snodgras/Fractionation/panand-annotations/.
```

Information about the methods for creating the annotations from Arun (6/14/22):
super brief prediction methods:
direct evidence from RNAseq assembled transcripts - finalized using mikado (basically fragmented transcripts to full length genes),
BRAKER for ab initio
for zea spp, sorghum as homology model donor
merged them using weights - higher weight to homology. When there is conflict - it resolves them using homology as a reference.
I thought you should know this, if you start noticing some bias.
3:21
This is very much like NAM, but they used full length transcpts to polish, which we don't have it here.

Also these gene ids are not comparable across genomes. ie: Zmay00001 identifies different genes across genomes, not the same gene

###1. Running Orthofinder
Make sure to use the `-og` flag
Need these input files

```
fasta of all protein coding gene's AA sequences(a proteome), uncompressed (not gzipped)
```

Recommendation from Orthofinder tutorial `https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html`


The files from Ensembl will contain many transcripts per gene. 
If we ran OrthoFinder on these raw files it would take 10x longer than necessary and could lower the accuracy. 
We’ll use a script provided with OrthoFinder to extract just the longest transcript variant per gene and run OrthoFinder on these files:

`for f in *fa ; do python ~/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done`

Shortening the filename is also a good ideas as it keeps the results tidy as the filenames are used to refer to the species, e.g. I shorten to Homo_sapiens.fa.


Need these output files
```
Blast*.txt
Orthogroups.tsv (Orthogroups.csv works to, be sure to specify this file using the -og flag)
SequenceIDs.txt
SpeciesIDs.txt
```

0. Reformat gff to primary transcript only and update gene IDs

The current (6/17/2022) gff files do not have genome specific gene IDs, just Zmays00001 or gene1 or mRNA1. 
They also have all isoforms for each gene. 
I will parse it down to the first (assumed longest) transcript for each gene, which is hopefully the most representative model for capturing exons.
I will also alter the gene IDs so that they are unique to each genome and thus won't cause downstream errors. 

`05.0.pSONIC.reformatGFF.sh`
This script will grep out all the primary transcripts and child features of the mRNA.
To get only the mRNA features, use awk '{print "ID="$2";Parent="$1}' on the `primary-transcript-ids.txt` file to get searchable strings
```
gff=$1
genome=$2

awk '$3=="mRNA"' ${gff} |cut -f 9 |cut -f 1-2 -d ";" |sed 's/;Parent=/\t/g' |sed 's/ID=_//g' |awk '{print $2"\t"$1}' | awk '!seen[$1]++' > ${genome}-primary-transcript-ids.txt
#awk '{print "ID="$2}' ${genome}-primary-transcript-ids.txt > IDsearchStrings.txt
cut -f 2  ${genome}-primary-transcript-ids.txt > IDsearchStrings.txt
echo "Beginning string search..."
grep -wF -f IDsearchStrings.txt ${gff} > ${genome}-primary-transcript.gff
echo "Finished string search. Deleting ID strings used for searching so there aren't mix-ups later..."
rm IDsearchStrings.txt

```
The slurm script (only takes a couple seconds to run)
```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

bash scripts/05.0.pSONIC.reformatGFF.sh panand-annotations/Zea-mays-ssp-parviglumis-TIL01_BIND.v1.gff3 Zv-TIL01
sed -i 's/=mRNA/=Zv-TIL01_mRNA/g' Zv-TIL01-primary-transcript.gff

bash scripts/05.0.pSONIC.reformatGFF.sh panand-annotations/Zea-mays-ssp-mexicana-TIL18_BIND.v1.gff3 Zx-TIL18
sed -i 's/=mRNA/=Zx-TIL18_mRNA/g' Zx-TIL18-primary-transcript.gff

```
Other ideas Arun suggested 
```
seqtk subseq full-pep.fa primary-transcript-ids.txt > primary-pep.fa

ml bioawk
bioawk -c fastx '{print ">ID_"$name"\n"$seq}' input.pep | fold > renamed_input.pep 

```

1. Convert gff to a proteome fasta `05.1.pSONIC.orthofinderInputs.sh`

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

Slurm script

```
#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=00:30:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 1 processor core(s) per node 
#SBATCH --mail-user=snodgras@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

grep scaf Zv-TIL01-primary-transcript.gff | awk -v OFS='\t' '$3 == "mRNA" {print $1}'| sort | uniq -c | tr ' ' '\t' | awk '$1 >= 5 {print $2}' - > Zv-TIL01.scaf.ids
grep -Fw -f Zv-TIL01.scaf.ids Zv-TIL01-primary-transcript.gff > Zv-TIL01-primary-5andup.gff
grep chr Zv-TIL01-primary-transcript.gff >> Zv-TIL01-primary-5andup.gff

grep scaf Zx-TIL18-primary-transcript.gff | awk -v OFS='\t' '$3 == "mRNA" {print $1}'| sort | uniq -c | tr ' ' '\t' | awk '$1 >= 5 {print $2}' - > Zx-TIL18.scaf.ids
grep -Fw -f Zx-TIL18.scaf.ids Zx-TIL18-primary-transcript.gff > Zx-TIL18-primary-5andup.gff
grep chr Zx-TIL18-primary-transcript.gff  >> Zx-TIL18-primary-5andup.gff

bash scripts/05.1.pSONIC.orthofinderInputs.sh panand-assemblies/Zv-TIL01-Reference-PanAnd-2.0.fasta Zv-TIL01-primary-5andup.gff Zv-TIL01
bash scripts/05.1.pSONIC.orthofinderInputs.sh panand-assemblies/Zx-TIL18-Reference-PanAnd-2.0.fasta Zx-TIL18-primary-5andup.gff Zx-TIL18
```

3. Run orthofinder
move all the fastas to `cds-fastas/` or `pep-fastas`

`slurm_05.2.pSONIC.runOrthofinder.sh` _Took ~13 minutes to run on 2 genomes_
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
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

ml gcc/10.2.0-zuvaafu
ml orthofinder/2.5.2-py3-z2m2odr
ml diamond

bioawk -c fastx '{gsub(/\./,"*",$seq); print ">"$name"\n"$seq}' input.fasta 

orthofinder -og -t $SLURM_JOB_CPUS_PER_NODE -f pep-fastas/ -S diamond
```

Notes on what to do once its finished
```
You should move these files into the same directory, and concatentate all of the Blast files into a single file named <PREFIX>.blast 
(where <PREFIX> is any name, but must be the same <PREFIX> used throughout).
```
Doing the above
```
mkdir TIL01-TIL18_pSONIC_test
cp pep-fastas/OrthoFinder/Results_Jun17/WorkingDirectory/SpeciesIDs.txt TIL01-TIL18_pSONIC_test/.
cp pep-fastas/OrthoFinder/Results_Jun17/WorkingDirectory/SequenceIDs.txt TIL01-TIL18_pSONIC_test/.
cp pep-fastas/OrthoFinder/Results_Jun17/Orthogroups/Orthogroups.tsv TIL01-TIL18_pSONIC_test/.
for i in pep-fastas/OrthoFinder/Results_Jun17/WorkingDirectory/Blast* ; do cp $i TIL01-TIL18_pSONIC_test/. ; done

cd TIL01-TIL18_pSONIC_test/
gunzip Blast*.gz
cat Blast*.txt > TIL01-TIL18.blast
```

My <PREFIX> will be TIL01-TIL18 for this test run

###2. Prepare files for MCScanX
install pSONIC
```
git clone https://github.com/conJUSTover/pSONIC
```

(recommended conda environment)
```
ml miniconda3
conda create -n pSONICenv python=3.7
source activate pSONICenv
pip install igraph 
```
Format gff files of all species into a "gff" file for MCScanX/pSONIC
* Need a two letter species code in front of each chromosome 
* Need the start and stop position of each gene
* A gene ID that is a unique  name for each gene in the fasta files

Want to also filter out any scaffolds that have less than 5 genes on it (MCScanX has 5 gene minimum and will slow down with a bunch of extra scaffolds)
Scaffolds just need to be numbered so that they start after the last chromosome.  

```
Sp##	GeneID	Start_POS	End_POS
```
VA = Zea mays spp parviglumis TIL01
XA  = Zea mays spp mexicana TIL18

`05.3.makeCombinedGFF.sh`
```
#gets them into the basic order
awk -v OFS='\t' '$3=="mRNA" {split($9,a,/;/); print $1, a[1], $4, $5}' Zv-TIL01-primary-transcript.gff > TIL01.temp.gff  
#removes the ID tag from the ID column
sed -i 's/ID=//g' TIL01.temp.gff
#Finds all the scaffolds that have more than 5 genes on them
#grep scaf TIL01.temp.gff | sort | cut -f 1 | uniq -c | tr ' ' '\t' | awk  '$1 >= 5 {print $2}' - > TIL01.scaf.ids
#adds all the chromosome positions to a different file
#grep chr TIL01.temp.gff >> TIL01.temp2.gff
#adds all the scaffolds with more than 5 genes to the different file
#grep -Fw -f TIL01.scaf.ids TIL01.temp.gff >> TIL01.temp2.gff 

awk -v OFS='\t' '$3=="mRNA" {split($9,a,/;/); print $1, a[1], $4, $5}' Zx-TIL18-primary-transcript.gff > TIL18.temp.gff  
sed -i 's/ID=//g' TIL18.temp.gff
#grep scaf TIL18.temp.gff | sort | cut -f 1 | uniq -c | tr ' ' '\t' | awk  '$1 >= 5 {print $2}' - > TIL18.scaf.ids
#grep chr TIL18.temp.gff >> TIL18.temp2.gff
#grep -Fw -f TIL18.scaf.ids TIL18.temp.gff >> TIL18.temp2.gff 

#remove the 'chr' and 'scaf_' characters, change the scaffold numbers
grep -w scaf_[1-9] *temp.gff
grep -w scaf_10 *temp.gff
#if the above doesn't return any lines, then there's no need to change the scaffold numbers

sed -i 's/chr/VA/g' TIL01.temp2.gff
sed -i 's/scaf_/VA/g' TIL01.temp2.gff
sed -i 's/chr/XA/g' TIL18.temp2.gff
sed -i 's/scaf_/XA/g' TIL18.temp2.gff

cat TIL01.temp2.gff TIL18.temp2.gff >> TIL01-TIL18.combined.gff
rm *temp* *scaf.ids
```

run the command to reformat the gff file
Should be named anything other than <PREFIX>.gff
`slurm_05.4.pSONIC.translategff.sh`

```
ml miniconda3
source activate pSONICenv
#python pSONIC.py <PREFIX> translate_gff -gff your_gff_file.gff

mv TIL01-TIL18.combined.gff TIL01-TIL18_pSONIC_test/.
cd TIL01-TIL18_pSONIC_test
python ../pSONIC/pSONIC.py TIL01-TIL18 translate_gff -gff TIL01-TIL18.combined.gff 

```

###3. Run MCScanX
install mcscanx
```
ml java (or jdk) #I can't remember which one was the one that ran make
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanX
make

#chmod +x MCScanX #Just need to make sure MCScanX is executable
```

`slurm_05.5.pSONIC.runMCScanX.sh`
```
ml miniconda3
source activate pSONICenv
cd TIL01-TIL18_pSONIC_test
../MCScanX/MCScanX -b 2 TIL01-TIL18
```

###4. Run pSONIC

`slurm_05.6.pSONIC.runpSONIC.sh`
```
ml miniconda3
source activate pSONICenv
cd TIL01-TIL18_pSONIC_test

python3 ../pSONIC/pSONIC.py TIL01-TIL18
```

##Running GENESPACE
"https://github.com/jtlovell/GENESPACE"

```
library(devtools)
usethis::create_github_token() #get new token
usethis::edit_r_environ(github_pat("R:GITHUB_PAT")) #copy token here

BiocManager::install(c("Biostrings", "rtracklayer"))
devtools::install_github("jtlovell/GENESPACE", upgrade = F)
#Do this in the shell (have load miniconda)
conda create -n orthofinder
source activate orthofinder
conda install -c bioconda orthofinder 

library(GENESPACE)
runwd <- file.path("/work/LAS/mhufford-lab/snodgras/Fractionation/genespace_test") #Where run should happen, empty or non-existant directory
```
Make sure that my datafiles are in the `runwd` in the following hierarchy
`[runWD]/rawGenomes/[speciesID]/[versionID]/annotation/[annotation gff and peptide fastas]`

```
#initialize GENESPACE run
gids <- c("Zv","Zx") #string list of genome names
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = c("TIL01","TIL18"), 
  ploidy = rep(1,2),
  wd = runwd, 
  gffString = "gff", 
  pepString = "pep",
  path2orthofinder = "~/.conda/envs/orthofinder/bin/orthofinder", #path to Orthofinder
  path2mcscanx = "/work/LAS/mhufford-lab/snodgras/Fractionation/MCScanX/", #path to MCScanX
  rawGenomeDir = file.path(runwd, "rawGenomes"))
  
#Format the raw annotations
#This may be tricky with my non ncbi annotations
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "gene", 
  gffIdColumn = "locus",
  gffStripText = "locus=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "locus=")
  
#Run Orthofinder
gpar <- run_orthofinder(gsParam = gpar)

#Run synteny search
gpar <- synteny(gsParam = gpar)

#Make multi-species riparian plot
ripdat <- plot_riparianHits(gpar)

#Build pangenome annotation
pg <- pangenome(gpar)
```

Current Trouble Shoot: 
```
> parse_annotations(gsParam = gpar, gffEntryType = "mRNA",gffIdColumn = "locus",gffStripText = "locus=", headerEntryIndex = 1,headerSep = " ",headerStripText = "locus=")
Parsing annotation files ...
	Zv ... 
		Importing gff ... found 37835 gff entires, and 37835 mRNA entries

**FOUND DUPLICATE GFF ENTRIES. THERE COULD BE A PROBLEM*** 
subsetting to longest model for each gene
		Importing fasta ... found 37835 fasta entires
		0 gff-peptide matches
Error in parse_annotations(gsParam = gpar, gffEntryType = "mRNA", gffIdColumn = "locus",  : 
  PARSING FAILED - run with troubleshoot = TRUE to see where the issues might be
  
Error in parse_annotations(gsParam = gpar, gffEntryType = "mRNA", gffIdColumn = "mRNA",  : 
  PARSING FAILED - try different parameters
```