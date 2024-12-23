# FRACTIONATION PIPELINE

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

###Running the full set of genomes
* Include all NAM genomes
* Include all teosinte genomes

1. Make the directory
_fullrun.1 is for defense, fullrun.2 includes TdKS_
```
mkdir fullrun.2
cd fullrun.2
mkdir genomeRepo
for i in TdKS TdFL ZdGigi ZdMomo ZhRIMHU001 ZlRIL003 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 Sb313 Av ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ; do mkdir genomeRepo/${i} ;done
```
There should be 38 genome directories in the `/genomeRepo`

2. Make the primary peptide fastas
Run `00.GENESPACE.primarytranscriptpep.sh` on genomes from scratch using the following slurm script:
```
cd fullrun.2/genomeRepo

#for pan-and genomes
for i in Td-KS Td-FL Zd-Gigi Zd-Momo Zh-RIMHU001 #Zl-RIL003 
	Zn-PI615697 Zv-TIL01 Zv-TIL11 Zx-TIL18 Zx-TIL25 Av; do
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/panand_annotations.v2/${i}*.gff3 ] 
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/panand_annotations.v2/${i}*.gff3.gz 
	fi
	if [ ! -f /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ]
	then gunzip /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta.gz
	fi
	bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/00.GENESPACE.primarytranscriptpep.sh /work/LAS/mhufford-lab/snodgras/Fractionation/panand_annotations.v2/${i}*.gff3 /work/LAS/mhufford-lab/snodgras/Fractionation/assemblies_final/${i}*.fasta ${i//\-}/${i//\-};
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
cd fullrun.2/genomeRepo/
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
Has to be done in an R environment

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

`slurm_02.minimapAnchors.sh`

## Filter Sb313 genes to be those that aren't repeats and found in Av

using Maggie's orthologs from NAM paper `Sbicolor_313_v3.1_non-tandem_gene_models_from_exons_primary_notandems_cshl_clusters4.txt`
and Arun's filtered orthogroups `sbicol_avirgi.tsv`

```
tail -n +2 Sbicolor_313_v3.1_non-tandem_gene_models_from_exons_primary_notandems_cshl_clusters4.txt | grep -w -f - sbicol_avirgi.tsv > test.txt


#filter CDS bed with Maggie's no tandem set
tail -n +2 ../code/Sbicolor_313_v3.1_non-tandem_gene_models_from_exons_primary_notandems_cshl_clusters4.txt | grep -w -f - Sb313.CDS.bed > notandems.Sb313.CDS.bed

#count the number of CDS (exons) per gene
##split the ID string into parent and child ID

cut -f 4 notandems.Sb313.CDS.bed | sed -e 's/;/\t/g' |cut -f 2 | sort | uniq -c | awk '{print $2"\t"$1}' - > Sb313.exonCntsPerGene.txt

#repeat with Av gff
awk -v OFS='\t' '$3 == "exon" && $7 == "+" {print $1,$4-1,$5,$9,$7}; $3 == "exon" && $7 == "-" {print $1,$4,$5+1,$9,$7}' Av-Kellogg1287_8-REFERENCE-PanAnd-1.0_final.gff3 > Av.exon.bed

cut -f 4 Av.exon.bed |sed -e 's/;/\t/g' | cut -f 1 | sort |uniq -c | awk '{print $2"\t"$1}' - > Av.exonCntsPerGene.txt

grep -v "," sbicol_avirgi.tsv | awk -v OFS='\t' 'NF == 3 {print $0}' - > one-one_sbicol_avirgi.tsv 

#Use the OGs from Arun to create a gene ID key by adding OG as a column to the exonCnts
while read line ; do 
	ID=$(echo $line | sed 's/Parent=//g' - |cut -f 1 -d " ")
	OG=$(grep -w $ID sbicol_avirgi.tsv | cut -f 1)
	echo $line $OG >> Av.exonCntsPerGene.OG.txt
done < Av.exonCntsPerGene.txt

while read line ; do 
	ID=$(echo $line | sed 's/Parent=//g' - | cut -f 1-3 -d ".")
	OG=$(grep -w $ID sbicol_avirgi.tsv | cut -f 1)
	echo $line $OG >> Sb313.exonCntsPerGene.OG.txt
done < Sb313.exonCntsPerGene.txt

while read line ; do
        ID=$(echo $line | sed 's/Parent=//g' - |cut -f 1 -d " ")
        OG=$(grep -w $ID one-one_sbicol_avirgi.tsv | cut -f 1)
        echo $line $OG >> Av.exonCntsPerGene.one-one.OG.txt
done < Av.exonCntsPerGene.txt

while read line ; do
        ID=$(echo $line | sed 's/Parent=//g' - | cut -f 1-3 -d ".")
        OG=$(grep -w $ID one-one_sbicol_avirgi.tsv | cut -f 1)
        echo $line $OG >> Sb313.exonCntsPerGene.one-one.OG.txt
done < Sb313.exonCntsPerGene.txt

#compare the number of exons per gene ID in Sb313 vs Av
sed -i "s/sbicol_avirgi.tsv://g" Av.exonCntsPerGene.OG.txt
sed -i "s/sbicol_avirgi.tsv://g" Av.exonCntsPerGene.one-one.OG.txt


join -1 3 -2 3 -e '-' -a1 -a2 <(sort -k 3 Av.exonCntsPerGene.OG.txt) <(sort -k 3 Sb313.exonCntsPerGene.OG.txt) > joined.exonCntsPerGene.OG.txt

join -1 3 -2 3 -e '-'  -a1 -a2 <(sort -k 3 Av.exonCntsPerGene.one-one.OG.txt) <(sort -k 3 Sb313.exonCntsPerGene.one-one.OG.txt) > joined.exonCntsPerGene.one-one.OG.txt

grep OG joined.exonCntsPerGene.one-one.OG.txt > joined.exonCntsPerGene.one-one.nonulls.OG.txt
awk -v OFS='\t' '$3 == $5 {print $0}' joined.exonCntsPerGene.one-one.nonulls.OG.txt > same.exonCntsPerGene.txt

wc -l same.exonCntsPerGene.txt 
9491 
awk 'NF == 5 {print}' joined.exonCntsPerGene.one-one.nonulls.OG.txt | wc -l 
15690

```
Maggie also created lift over files of each Sb annotation and Av (reciprocal liftovers)
```
unzip Av_sorghum_lifted_gff_files_and_CDS_counts.zip 
unzip Av-Kellogg1287_8-REFERENCE-PanAnd-1.0_noexon_lifted_to_Sbicolor_313_v3.1.gene_synteny_status_sorted.gff3.zip 

cd Av_sorghum_lifted_gff_files_and_CDS_counts/counts_of_CDSs_that_match_in_length_and_location/

awk -v OFS='\t' '$4 == $5 {print $0}' Av-Kellogg1287_8-REFERENCE-PanAnd-1.0_noexon_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt | wc -l 
  21081

wc -l Av-Kellogg1287_8-REFERENCE-PanAnd-1.0_noexon_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt 
   79716 Av-Kellogg1287_8-REFERENCE-PanAnd-1.0_noexon_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt
   
awk -v OFS='\t' '$4 == $5 {print $0}' Av-Kellogg1287_8-REFERENCE-PanAnd-1.0_noexon_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt | awk -v OFS='\t' '{tmp = $5 / $2 }; {print $1,$2,$3,$4,tmp}' - > Av_Sb313_CDS-counts.allCDSsamelength.percentage.txt
#This give the counts for a histogram of the coverage of the exact matched Av exons to the Sb exons
cut -f 5 Av_Sb313_CDS-counts.allCDSsamelength.percentage.txt | sort | uniq -c 
#notably, of the 21K, ~18K have exact match (Same number of exons, exact same lengths) (exact number is 17933)

awk -v OFS="\t" '$5 == 1 {print $0 }' Av_Sb313_CDS-counts.allCDSsamelength.percentage.txt | cut -f 1 | sort | uniq | wc -l  
   17906 #number of exact matches that are unique Sb313 gene IDs

awk -v OFS="\t" '$5 == 1 {print $0 }' Av_Sb313_CDS-counts.allCDSsamelength.percentage.txt > Av_Sb313_ExactMatchesOnly.txt

grep -f ../../Sbicolor_313_v3.1_non-tandem_gene_models_from_exons_primary_notandems_cshl_clusters4.txt Av_Sb313_ExactMatchesOnly.txt > Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt
wc -l Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt 
   15111 Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt

cut -f 1 Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt | cut -f 1-2 -d "." | sort | uniq | wc -l 
   12169 #Just uniq gene models?
   
#How to deal with isoforms:
#keep the one with the most exons
#if a tie, pick one at random
#This will be easiest to do in R
   
#Let's repeat with the other sorghum annotations:
#SbicolorBTx642_564_v1.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt
#SbicolorRTx430_552_v2.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt
#SbicolorSC187_694_v1.1.gene_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt

#how many where there is an exact match between the number of lifted exons and the number of exons with same length
for i in SbicolorBTx642_564_v1.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt SbicolorRTx430_552_v2.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt SbicolorSC187_694_v1.1.gene_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt ; do
awk -v OFS='\t' '$4 == $5 {print $0}' $i | wc -l 
done

   36287
   36098
   36449

#how many lifted over gene models total
for i in SbicolorBTx642_564_v1.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt SbicolorRTx430_552_v2.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt SbicolorSC187_694_v1.1.gene_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt ; do
wc -l $i
done
   77059 SbicolorBTx642_564_v1.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt
   72157 SbicolorRTx430_552_v2.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt
   75738 SbicolorSC187_694_v1.1.gene_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt

#make percentage column for number of lifted over 100% exact match exons/# of Sb313 exons in model
for i in SbicolorBTx642_564_v1.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt SbicolorRTx430_552_v2.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt SbicolorSC187_694_v1.1.gene_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt ; do
awk -v OFS='\t' '$4 == $5 {print $0}' $i | awk -v OFS='\t' '{tmp = $5 / $2 }; {print $1,$2,$3,$4,tmp}' - > ${i%_v1.1.gene*txt}_Sb313.allCDSsamelength.percentage.txt
done
mv SbicolorRTx430_552_v2.1.gene_Chr_vs_Sbicolor_313_v3.1.gene_CDS-counts.txt_Sb313.allCDSsamelength.percentage.txt SbicolorRTx430_552_Sb313.allCDSsamelength.percentage.txt

#This gets the number for exact match instead of creating all histogram values like the uniq -c command did for Av
for i in Sbicolor*percentage.txt ; do 
cut -f 5 $i | awk -v OFS='\t' '$1 == 1 {print $0}' - | wc -l
done
   35558 #97.99% of original exact match
   35529 #98.42% of original exact match
   35909 #98.52% of original exact match

#This removes multimapping to the same gene model
for i in Sbicolor*percentage.txt ; do awk -v OFS='\t' '$5 == 1 {print $0}' $i | cut -f 1 | sort | uniq | wc -l; done
   28186
   28375
   26655


```
Collapse exons 
```
bedtools sort -i test.Sb313.CDS.bed | bedtools merge -c 4 -o distinct -i - > merged.test.Sb313.CDS.bed

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

03.runAnchorWave.ZnZd.sh
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
	-R 4 \ 
	-Q 1 \
	-o AnchorWave_output/Sb313_${GenomeName}_anchorwave.4to1.maf \ 
	-t 5 > AnchorWave_logs/Sb313_${GenomeName}_anchorwave.4to1.log 2>&1 
```
```
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a TdFL >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Av.sh assemblies_final/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0 Av  >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.ZnZd.sh assemblies_final/Zn-PI615697-REFERENCE-PanAnd-1.0 ZnPI615697 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zv-TIL01-REFERENCE-PanAnd-1.0 ZvTIL01 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.ZnZd.sh assemblies_final/Zd-Gigi-REFERENCE-PanAnd-1.0 ZdGigi >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zv-TIL11-REFERENCE-PanAnd-1.0 ZvTIL11 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.ZnZd.sh assemblies_final/Zd-Momo-REFERENCE-PanAnd-1.0 ZdMomo >> scripts/AnchorWave.cmds.txt    
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zx-TIL18-REFERENCE-PanAnd-1.0 ZxTIL18 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zh-RIMHU001-REFERENCE-PanAnd-1.0 ZhRIMHU001 >> scripts/AnchorWave.cmds.txt      
echo bash scripts/03.runAnchorWave.Tripsacinae.sh assemblies_final/Zx-TIL25-REFERENCE-PanAnd-1.0 ZxTIL25 >> scripts/AnchorWave.cmds.txt

echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-B73 ZmB73 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-B97 ZmB97 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML103 ZmCML103 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML228 ZmCML228 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML247 ZmCML247 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML277 ZmCML277 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML322 ZmCML322 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML333 ZmCML333 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML52 ZmCML52 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-CML69 ZmCML69 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-HP301 ZmHP301 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-IL14H ZmIL14H >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Ki11 ZmKi11 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Ki3 ZmKi3 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Ky21 ZmKy21 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-M162W ZmM162W >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-M37W ZmM37W >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Mo18W ZmMo18W >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-MS71 ZmMS71 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-NC350 ZmNC350 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-NC358 ZmNC358 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Oh43 ZmOh43 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Oh7b ZmOh7b >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-P39 ZmP39 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Tx303 ZmTx303 >> scripts/AnchorWave.cmds.txt
echo bash scripts/03.runAnchorWave.Tripsacinae.sh NAM-assemblies/Zm-Tzi8 ZmTzi8 >> scripts/AnchorWave.cmds.txt

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
```

## Identifying gap regions from the assemblies

```
#Testing on Zd-Gigi first
cd assemblies_final/
module load seqtk
seqtk cutN -gp10000000 -n1 Zd-Gigi-REFERENCE-PanAnd-1.0.fasta > Zd-Gigi.ngaps.bed

grep "chr9" ZdGigi_chr9_Chr10_sorted_filtered.maf | cut -f 1-7 -d " " | awk -v OFS='\t' '$4 > $3 {print $2,$3,$4}; $4 < $3 {print $2,$4,$3}' - > Zd-Gigi-maftemp.bed

ml bedtools2
bedtools intersect -a Zd-Gigi-maftemp.bed -b ../assemblies_final/Zd-Gigi.ngaps.txt 

for i in Av-Kellogg1287_8 Td-FL Zd-Momo Zh-RIMHU001 Zl-RIL003 Zn-PI615697 Zv-TIL01 Zv-TIL11 Zx-TIL18 Zx-TIL25 ; do
seqtk cutN -gp10000000 -n1 ${i}*.fasta > ${i}.ngaps.bed ; done

cd ../NAM-assemblies/
for i in *.fasta ; do
seqtk cutN -gp10000000 -n1 ${i} > ${i%.fasta}.ngaps.bed ; done 

```
##Testing to see if the gap regions are causing issues with the deletion calls
```
mkdir unzipped_gvcf
gunzip trimmed_Sb313_ZdGigi_Chr10.bin1.gvcf.gz 
gunzip trimmed_Sb313_ZdGigi_Chr10.bin2.gvcf.gz 
mv trimmed*gvcf unzipped_gvcf/
cd unzipped_gvcf
grep -v "^#" trimmed_Sb313_ZdGigi_Chr10.bin1.gvcf | awk -v OFS='\t' '{print $8,$10}' - > ZdGigi.bin1.genotypes.txt
cut -f 2 ZdGigi.bin1.genotypes.txt | cut -f 1 -d ":" | sort | uniq -c
grep -v "^#" trimmed_Sb313_ZdGigi_Chr10.bin2.gvcf | awk -v OFS='\t' '{print $8,$10}' - > ZdGigi.bin2.genotypes.txt
cut -f 2 ZdGigi.bin2.genotypes.txt | cut -f 1 -d ":" | sort | uniq -c

cut -f 1 ZdGigi.bin1.genotypes.txt | cut -f 1 -d ";" - | sed "s/ASM_Chr=/chr/g" > bin1.chr.temp
cut -f 1 ZdGigi.bin1.genotypes.txt | cut -f 2 -d ";" - | sed "s/ASM_End=//g"> bin1.end.temp
cut -f 1 ZdGigi.bin1.genotypes.txt | cut -f 3 -d ";" - | sed "s/ASM_Start=//g" > bin1.start.temp

cut -f 1 ZdGigi.bin1.genotypes.txt | paste bin1.chr.temp bin1.start.temp bin1.end.temp - > bin1.bed.temp
awk -v OFS='\t' '$2 < $3 {print $0}; $2 > $3 {print $1, $3, $2, $4}' bin1.bed.temp > bin1.bed

ml bedtools2
bedtools intersect -wb -a ../../assemblies_final/Zd-Gigi.ngaps.txt -b bin1.bed | cut -f 7 - > bin1.gaps.to.check.txt
grep -f bin1.gaps.to.check.txt trimmed_Sb313_ZdGigi_Chr10.bin1.gvcf > bin1.gap.vcf.lines

while read line ; do echo $line | cut -f 5 -d " "| wc -c ; done <bin1.gap.vcf.lines 
#This counts the number of characters for the ALT allele
while read line ; do echo $line | cut -f 4 -d " " |wc -c ; done <bin1.gap.vcf.lines 
#This counts the number of characters for the REF allele
####NOTE THAT ALL THE ONES THAT INTERSECT FOR BIN1 ARE INSERTIONS

# repeat for bin2
cut -f 1 ZdGigi.bin2.genotypes.txt | cut -f 1 -d ";" - | sed "s/ASM_Chr=/chr/g" > bin2.chr.temp
cut -f 1 ZdGigi.bin2.genotypes.txt | cut -f 2 -d ";" - | sed "s/ASM_End=//g"> bin2.end.temp
cut -f 1 ZdGigi.bin2.genotypes.txt | cut -f 3 -d ";" - | sed "s/ASM_Start=//g" > bin2.start.temp

cut -f 1 ZdGigi.bin2.genotypes.txt | paste bin2.chr.temp bin2.start.temp bin2.end.temp - > bin2.bed.temp
awk -v OFS='\t' '$2 < $3 {print $0}; $2 > $3 {print $1, $3, $2, $4}' bin2.bed.temp > bin2.bed

bedtools intersect -wb -a ../../assemblies_final/Zd-Gigi.ngaps.txt -b bin2.bed | cut -f 7 - > bin2.gaps.to.check.txt
grep -f bin2.gaps.to.check.txt trimmed_Sb313_ZdGigi_Chr10.bin2.gvcf > bin2.gap.vcf.lines

while read line ; do echo $line | cut -f 5 -d " "| wc -c ; done <bin2.gap.vcf.lines 
while read line ; do echo $line | cut -f 4 -d " " |wc -c ; done <bin2.gap.vcf.lines 

##Turn it into a bash script (checkgaps.sh)
ngap=$1 #../../assemblies_final/Zd-Gigi.ngaps.txt
gvcf=$2 #trimmed_Sb313_ZdGigi_Chr10.bin1.gvcf.gz
out=$3 #ZdGigi.bin1

if [[ ! -f ${gvcf%.gz} ]] ; then gunzip ${gvcf} ; fi

grep -v "^#" ${gvcf%.gz} |  awk -v OFS='\t' '{print $8,$10}' - > ${out}.genotypes.txt
cut -f 1 ${out}.genotypes.txt | cut -f 1 -d ";" - | sed "s/ASM_Chr=/chr/g" > ${out}.chr.temp
cut -f 1 ${out}.genotypes.txt | cut -f 2 -d ";" - | sed "s/ASM_End=//g"> ${out}.end.temp
cut -f 1 ${out}.genotypes.txt | cut -f 3 -d ";" - | sed "s/ASM_Start=//g" > ${out}.start.temp

cut -f 1 ${out}.genotypes.txt | paste ${out}.chr.temp ${out}.start.temp ${out}.end.temp - > ${out}.bed.temp
awk -v OFS='\t' '$2 < $3 {print $0}; $2 > $3 {print $1, $3, $2, $4}' ${out}.bed.temp > ${out}.bed

ml bedtools2
bedtools intersect -wb -a ${ngap} -b ${out}.bed | cut -f 7 - > ${out}.gaps.to.check.txt
grep -f ${out}.gaps.to.check.txt ${gvcf%.gz} > ${out}.gap.vcf.lines

echo The number of characters for ${out} ALT alleles with gaps: 
while read line ; do echo $line | cut -f 5 -d " "| wc -c ; done <${out}.gap.vcf.lines
echo The number of characters for ${out} REF alleles with gaps: 
while read line ; do echo $line | cut -f 4 -d " " |wc -c ; done <${out}.gap.vcf.lines 


bash checkgaps.sh ../../assemblies_final/Td-FL.ngaps.bed ../trimmed_Sb313_TdFL_Chr10.bin1.gvcf.gz TdFL.bin1
bash checkgaps.sh ../../assemblies_final/Td-FL.ngaps.bed ../trimmed_Sb313_TdFL_Chr10.bin2.gvcf.gz TdFL.bin2

bash checkgaps.sh ../../assemblies_final/Zv-TIL01.ngaps.bed ../trimmed_Sb313_ZvTIL01_Chr10.bin1.gvcf.gz ZvTIL01.bin1
bash checkgaps.sh ../../assemblies_final/Zv-TIL01.ngaps.bed ../trimmed_Sb313_ZvTIL01_Chr10.bin2.gvcf.gz ZvTIL01.bin2

bash checkgaps.sh ../../NAM-assemblies/Zm-B73.ngaps.bed ../trimmed_Sb313_ZmB73_Chr10.bin1.gvcf.gz ZmB73.bin1
bash checkgaps.sh ../../NAM-assemblies/Zm-B73.ngaps.bed ../trimmed_Sb313_ZmB73_Chr10.bin2.gvcf.gz ZmB73.bin2
```
_It looks like any gap regions are always within large insertion calls where the ALT allele is much larger than the REF allele_

## Swap coordinates of MAF
Use script `https://github.com/baoxingsong/AnchorWave/blob/master/scripts/anchorwave-maf-swap.py`
```
genome=$1
cat Sb313_${genome}_anchorwave.maf | python /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/anchorwave-maf-swap.py > Sb313_${genome}_anchorwave.swap.maf

###SLURM
bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/03.5.swapMAF.sh TdFL

for i in *.maf ; do n=$(echo $i | cut -f 2 -d "_" ) ; echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/03.5.swapMAF.sh ${n} >> ../scripts/swap.cmds.txt ; done
#remove the first line so no swapping of Av (no need to split it)

#make sure there's 3 hours of time
for i in {1..35} ; do sbatch --dependency=afterok:4781250 ../scripts/swap.cmds_${i}.sub ;done

### Redoing Zd and Zn genomes with .4to1 + TdKS
#for i in *4to1.maf ; do do n=$(echo $i | cut -f 2 -d "_" ) ; echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/03.5.swapMAF.sh ${n}.4to1 >> ../scripts/swap.cmds.txt ; done
#Won't work because the 4to1 in file name will mess up the script
echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/03.5.swapMAF.sh TdKS >> ../scripts/swap.cmds.txt
```

#Split MAF (by parsing coordinates?)

*NOTE: There is no chr6 in the bin coordinates for TdFL*
At this point, manually curate the associations between tripsacinae chromosomes, subgenome identity, and sorghum chromosome using the calls of B73 in Schnable et al 2011

*kentutils*
`/ptmp/arnstrm/kentutils.sif`

see `split_and_filter_maf_files_for_Sam.txt`
First need to swap mafs before doing this
```mafsplit.1.sh
#!/bin/bash
ml singularity
genome=$1
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName swap_${genome}/ Sb313_${genome}_anchorwave.swap.maf

###
for i in *swap.maf ; do 
	n=$(echo $i | cut -f 2 -d "_")
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/mafsplit.1.sh $n >> ../scripts/mafsplit.cmds.txt 
done

cd ../scripts
python makeSLURM.py 1 mafsplit.cmds.txt
cd -
sbatch ../scripts/mafsplit.cmds_0.sub
for i in {1..35} ; do sbatch --dependency=afterok:4781359 ../scripts/mafsplit.cmds_${i}.sub ; done

singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName swap_ZnPI615697_4to1/ Sb313_ZnPI615697.4to1_anchorwave.swap.maf 
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName swap_ZdGigi_4to1/ Sb313_ZdGigi.4to1_anchorwave.swap.maf 
singularity exec --bind $PWD /ptmp/arnstrm/kentutils.sif mafSplit -byTarget dummy.bed -useFullSequenceName swap_ZdMomo_4to1/ Sb313_ZdMomo.4to1_anchorwave.swap.maf 

###

for f in */*.maf ;do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp"_"$fl"; done
rm *scaf*
rm *alt-ctg*

for f in swap_*4to1/*.maf swap_TdKS/*.maf ; do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp"_"$fl"; done
rm *scaf*
rm *ctg*
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

for f in TdKS_*/*.maf Z*4to1_chr*/*.maf ; do fp=$(dirname "$f"); fl=$(basename "$f"); mv "$fp/$fl" "$fp"_"$fl"; done
```

Manually:
```
mkdir tripsacinae-sb_split_mafs
mv *chr*_Chr*.maf tripsacinae-sb_split_mafs
cd tripsacinae-sb_split_mafs
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
mkdir test_mafs/
cat TdFL_chr13_Chr10_sorted_filtered.maf > test_mafs/TdFL_Chr10.bin2.maf
cat TdFL_chr10_Chr10_sorted_filtered.maf > test_mafs/TdFL_Chr10.bin1.maf
tail -n +2 TdFL_chr4_Chr10_sorted_filtered.maf >> test_mafs/TdFL_Chr10.bin1.maf

cat ZdGigi_chr5_Chr10_sorted_filtered.maf  > test_mafs/ZdGigi_Chr10.bin1.maf
tail -n +2 ZdGigi_chr9_Chr10_sorted_filtered.maf  >> test_mafs/ZdGigi_Chr10.bin1.maf
cat ZdGigi_chr6_Chr10_sorted_filtered.maf  > test_mafs/ZdGigi_Chr10.bin2.maf

cat ZvTIL01_chr5_Chr10_sorted_filtered.maf  > test_mafs/ZvTIL01_Chr10.bin1.maf
tail -n +2 ZvTIL01_chr9_Chr10_sorted_filtered.maf  >> test_mafs/ZvTIL01_Chr10.bin1.maf
cat ZvTIL01_chr6_Chr10_sorted_filtered.maf  > test_mafs/ZvTIL01_Chr10.bin2.maf

cat ZmB73_chr5_Chr10_sorted_filtered.maf  > test_mafs/ZmB73_Chr10.bin1.maf
tail -n +2 ZmB73_chr9_Chr10_sorted_filtered.maf  >> test_mafs/ZmB73_Chr10.bin1.maf
cat ZmB73_chr6_Chr10_sorted_filtered.maf  > test_mafs/ZmB73_Chr10.bin2.maf

wc -l *Chr10_sorted_filtered.maf #668
wc -l *Chr10.bin*maf #664 (4 less because removing header line of 1 maf file)
```
Merge split mafs into bins for whole data set by Sb313 chromosome
*Make sure the parsing files have a blank line at the end for this loop because otherwise, the concatenation won't work*
```
#make each bin file per query genome first (by Sb313 chromosome)

##generate unique list of query chromosome names per bin
#Using the results.v2 files
for i in zmB73vs*_results.tsv; do
	n=$(echo $i | sed 's/zmB73vs//g' | cut -f 1 -d "_")
	cut -f 1,7,9 $i | sort | uniq >> ${n}.bins.txt
done

##concatenate
#Do a while loop that uses to the two fields
slurm_03.1.concatMAFtoBins.sh
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 36
#SBATCH --mem=350GB
#SBATCH -t 3:00:00
#SBATCH -J featureCounts.
#SBATCH -o featureCounts.o%j
#SBATCH -e featureCounts.e%j
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

for g in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	while read -r query sb bin ; do 
		if [ -f "$g"_"$sb"."$bin".maf ] ; then 
			tail -n +2 "$g"_"$query"_"$sb"_sorted_filtered.maf >> "$g"_"$sb"."$bin".maf
			else
			cat "$g"_"$query"_"$sb"_sorted_filtered.maf > "$g"_"$sb"."$bin".maf
		fi
	done < /work/LAS/mhufford-lab/snodgras/Fractionation/parsing_coordinates/${g}.bins.txt
done

###
#problem, the bin included the return character into the name

for g in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
for i in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 ; do 
	mv ${g}_${i}.bin1*.maf ${g}_${i}.bin1.maf
	mv ${g}_${i}.bin2*.maf ${g}_${i}.bin2.maf
done
done


##double check that this works before launching on all
md5sum ZdGigi_Chr10.bin*maf
md5sum test_mafs/ZdGigi_Chr10.bin*maf

#Didn't make the TdFL and had an issue with some of the renaming of the Oh7b files
#ZmOh7b_Chr01.bin2 ZmOh7b_Chr07.bin1.maf

rm ZmOh7b_Chr01.bin2*maf ZmOh7b_Chr07.bin1*maf

cat ZmOh7b_chr5_Chr01_sorted_filtered.maf > ZmOh7b_Chr01.bin2.maf
tail -n +2 ZmOh7b_chr9_Chr01_sorted_filtered.maf >> ZmOh7b_Chr01.bin2.maf

cat ZmOh7b_chr1_Chr07_sorted_filtered.maf > ZmOh7b_Chr07.bin1.maf
tail -n +2 ZmOh7b_chr10_Chr07_sorted_filtered.maf >> ZmOh7b_Chr07.bin1.maf
tail -n +2 ZmOh7b_chr6_Chr07_sorted_filtered.maf >> ZmOh7b_Chr07.bin1.maf
tail -n +2 ZmOh7b_chr9_Chr07_sorted_filtered.maf >> ZmOh7b_Chr07.bin1.maf

```
Just remaking the Chr10 bin1 mafs due to issue with last line not being read by the loop above
```
for g in ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
cat "$g"_chr5_Chr10_sorted_filtered.maf > "$g"_Chr10.bin1.maf
tail -n +2 "$g"_chr9_Chr10_sorted_filtered.maf >> "$g"_Chr10.bin1.maf
echo Done with $g
done
```
For TdKS and the 4to1 redos of Zn and Zd
```
while read -r query sb bin ; do 
if [ -f TdKS_"$sb"."$bin".maf ] ; then 
			tail -n +2 TdKS_"$query"_"$sb"_sorted_filtered.maf >> TdKS_"$sb"."$bin".maf
			else
			cat TdKS_"$query"_"$sb"_sorted_filtered.maf > TdKS_"$sb"."$bin".maf
		fi
	done < /work/LAS/mhufford-lab/snodgras/Fractionation/parsing_coordinates/TdFL.bins.txt

for g in ZdGigi_4to1 ZdMomo_4to1 ZnPI615697_4to1 ; do 
while read -r query sb bin ; do 
		if [ -f "$g"_"$sb"."$bin".maf ] ; then 
			tail -n +2 "$g"_"$query"_"$sb"_sorted_filtered.maf >> "$g"_"$sb"."$bin".maf
			else
			cat "$g"_"$query"_"$sb"_sorted_filtered.maf > "$g"_"$sb"."$bin".maf
		fi
	done < /work/LAS/mhufford-lab/snodgras/Fractionation/parsing_coordinates/${g%_4to1}.bins.txt
done
```
note: there was no TdKS chr4 to Chr01 or Chr10

## Converting to GVCF and see if it works (use the `test_mafs`)
make sure that `--fillGaps false` is in the options for the plugin
```04.splitmaf2gvcf.sh
#!/bin/bash

REFfasta=$1 #Sorghum genome fasta
MAF=$2 #query name without extensions

##Convert maf to gvcf (this takes about 15 or so minutes):

#MAFToGVCFPlugin <options>
#-referenceFasta <Reference Fasta> : Input Reference Fasta (required)
#-mafFile <Maf File> : Input MAF file.  Please note that this needs to be a MAF file with 2 samples.  The first will 
be assumed to be the Reference and the second will be the assembly. (required)
#-sampleName <Sample Name> : Sample Name to write to the GVCF file as the genome header or ID (required)
#-gvcfOutput <Gvcf Output> : Output GVCF file name (required)
#-fillGaps <true | false> : When true, if the maf file does not fully cover the reference genome any gaps in coverage will be #filled in with reference blocks. This is necessary if the resulting GVCFs are to be combined. (Default: false)

module load singularity
singularity exec --bind $PWD /work/LAS/mhufford-lab/elena/.conda/phg_latest.sif /tassel-5-standalone/run_pipeline.pl \
        -Xmx300g \
        -MAFToGVCFPlugin \
        -referenceFasta ${REFfasta} \
        -mafFile ${MAF} \
        -sampleName Sb313_${MAF%.maf} \
        -gvcfOutput Sb313_${MAF%.maf}.gvcf \
        -fillGaps false > ${MAF%.maf}.gvcf.log 2>&1
        

##The command "-Xmx300g" demands that Java has enough RAM for the job to be run; in this case, 300g
```

```
for i in *bin*.maf ; do 
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/04.splitmaf2gvcf.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.fasta ${i} >> splitmaf2gvcf.cmds.txt
done
python ../../../scripts/makeSLURM.py 1 splitmaf2gvcf.cmds.txt

#test before running all
sbatch splitmaf2gvcf.cmds_0.sub 

for i in {1..7} ; do sbatch splitmaf2gvcf.cmds_${i}.sub ; done
```

For the full set of genomes
```
for i in *bin*.maf ; do 
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/04.splitmaf2gvcf.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.fasta ${i} >> splitmaf2gvcf.cmds.txt
done
python ../../scripts/makeSLURM.py 1 splitmaf2gvcf.cmds.txt

sbatch splitmaf2gvcf.cmds_0.sub 

for i in {1..699} ; do sbatch splitmaf2gvcf.cmds_${i}.sub ;done
```
For redoing just Chr10 bin1 mafs for Zea (not TdFl)
```
for i in Z*_Chr10.bin1.maf ; do 
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/04.splitmaf2gvcf.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.fasta ${i} >> Chr10redo2gvcf.cmds.txt
done
python ../../scripts/makeSLURM.py 1 Chr10redo2gvcf.cmds.txt
```
For TdKS and 4to1 Zn and Zd genomes
```
for i in TdKS*bin*.maf Z*4to1*bin*.maf ; do 
	echo bash /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/04.splitmaf2gvcf.sh /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.fasta ${i} >> TdKSandZnZd2gvcf.cmds.txt
done
python ../../scripts/makeSLURM.py 1 TdKSandZnZd2gvcf.cmds.txt
```

Run on the full set of genomes
_Same slurm script as above, just edit the time to be 24 hours (just in case)_

convert to vcf
Let's try keeping each bin separate:
(use the `test_mafs` first)

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

_Note: since we're using a limited reference set of homoeologs for the CDS_
I wrote the `ref_Sb313.cds` object from the R script `trial-fractionation-calling.R` to `/work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed`


To be able to load the vcf into R, must remove the header information
It will also speed things up to filter for deletions that cover ref exons

#####THIS IS THE PART THAT NEEDS TO BE REVISED
1. Genome wide measure of how much alignment is there between target genome and sorghum?
a. in general
b. by chromosome
c. how much of these alignments are indels? 

### Option 1: VCF metrics plugin from tassel 
This should give metrics of the GVCF, but I'm not entirely sure what those metrics are since I can't find documentation

_You have to move all the GVCFs into a single GVCF folder to work_

```
ml singularity
singularity exec --bind $PWD /ptmp/arnstrm/phg_latest.sif /tassel-5-standalone/run_pipeline.pl -VCFMetricsPlugin \
        	-vcfDir /work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/GVCF/ \
        	-outFile /work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/QC/GVCFmetrics.tsv \
        	-endPlugin
```


### Option 2: Look at read depth in a SAM made from the MAFs

```
#**Michelle’s easiest solution after thinking:
#** inspired by Baoxing’s supplement

## convert maf to sorted bam
maf-convert sam anchorwave.maf | sed 's/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.//g' | sed 's/Zea_mays.AGPv4.dna.toplevel.fa.//g' | samtools view -O BAM --reference Zea_mays.AGPv4.dna.toplevel.fa - | samtools sort - > anchorwave.bam


## count aligned sites (nucleotide-gap or nucleotide-nucleotide alignment)
samtools depth alignment.bam | wc -l

## count position match (nucleotide-nucleotide alignment, NOT counting whether they're snps!!)
samtools depth alignment.bam | awk ‘$3>0 {print $0}’ | wc -l


## count aligned sites in specific bed file
## if you don't do the wc -l, you'll get each site alignment stats
samtools depth alignment.bam -b all_reproducible_peaks_summits_merged.bed | wc -l

## count position match sites in specific bed file
samtools depth alignment.bam -b all_reproducible_peaks_summits_merged.bed | awk '$3>0{print $0}' | wc -l

## depth will be 2 for polyploids (-R2) - it may be useful to split "reads" aka alignment blocks into separate files to keep the phasing
## ask me more if this ends up being an issue!

```
Get the `maf-convert.py` from here: `https://gitlab.com/mcfrith/last/-/blob/main/bin/maf-convert`

How this translates to these files:
```
#Use the unswapped, because otherwise the SAM won't be read properly
#none of the strings to remove like the above so just funnel it samtools

ml samtools
python maf-convert.py sam AnchorWave_output/Sb313_TdFL_anchorwave.maf | samtools view -O BAM --reference Sb313.fasta - | samtools sort - >  AnchorWave_maf2bam/Sb313_TdFL_anchorwave.bam

### Expected length of the Sb genome:
ml bioawk
bioawk -c fastx '{ print $name, length($seq) }' < Sb313.fasta
Chr01	80884392
Chr02	77742459
Chr03	74386277
Chr04	68658214
Chr05	71854669
Chr06	61277060
Chr07	65505356
Chr08	62686529
Chr09	59416394
Chr10	61233695

bioawk -c fastx '{ print $name, length($seq) }' < Sb313.fasta | awk -F '\t' 'BEGIN{SUM=0}{ SUM+=$2 }END{print SUM}' - 
732152042
#for 2 complete copies its: 1464304084

## count aligned sites (nucleotide-gap or nucleotide-nucleotide alignment)
samtools depth AnchorWave_maf2bam/Sb313_TdFL_anchorwave.bam | wc -l
###What's the expectation here? 2x the length of the Sb genome?###

## count position match (nucleotide-nucleotide alignment, NOT counting whether they're snps!!)
samtools depth AnchorWave_maf2bam/Sb313_TdFL_anchorwave.bam | awk ‘$3>0 {print $0}’ | wc -l
###What's the expectation here?###

###Expected number of base pairs in the exons:
cat ref_Sb313.cds.bed | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
16390866

## count aligned sites in specific bed file
## if you don't do the wc -l, you'll get each site alignment stats
awk -v OFS='\t' '{print "Chr"$0}' ref_Sb313.cds.bed | samtools depth  AnchorWave_maf2bam/Sb313_TdFL_anchorwave.bam -b - | wc -l

## count position match sites in specific bed file
awk -v OFS='\t' '{print "Chr"$0}' ref_Sb313.cds.bed | samtools depth  AnchorWave_maf2bam/Sb313_TdFL_anchorwave.bam -b - | awk '$3>0{print $0}' | wc -l

## depth will be 2 for polyploids (-R2) - it may be useful to split "reads" aka alignment blocks into separate files to keep the phasing

```
Make it into a loop:
```AnchorWave_output/QC/slurm-maf2bam.sh
ml samtools
cd /work/LAS/mhufford-lab/snodgras/Fractionation/
for i in AnchorWave_output/Sb313_Z*anchorwave.maf ; do
        python maf-convert.py sam ${i} | samtools view -O BAM --reference Sb313.fasta - | samtools sort - >  AnchorWave_maf2bam/${i#Anchorwave_output/}.bam
done
```
Create a depth metrics file
```
ml samtools
echo Genome Sb_cnt_aligned_sites Sb_cnt_matched_sites Exon_cnt_aligned_sites Exon_cnt_matched_sites > depth.metrics.txt
for i in *.bam ; do 
	sbd=$(samtools depth $i | wc -l ) #sorghum depth of aligned regions
	sbm=$(samtools depth $i | awk '$3>0 {print $0}' | wc -l ) #sorghum depth of matched regions
	ed=$(awk -v OFS='\t' '{print "Chr"$0}' ../ref_Sb313.cds.bed | samtools depth $i -b - | wc -l ) #ref exon depth of aligned regions
	em=$(awk -v OFS='\t' '{print "Chr"$0}' ../ref_Sb313.cds.bed | samtools depth $i -b - | awk '$3>0{print $0}' | wc -l ) #ref exon depth of matched regions
	g=$(echo $i | cut -f 2 -d "_")
	echo $g $sbd $sbm $ed $em >> depth.metrics.txt
done
cat depth.metrics.txt | tr ' ' '\t' > depth.metrics.tsv

#Make depth files so I can look at average depth and not just number of positions aligning
for i in *.bam ; do 
samtools depth $i > ${i%.bam}.depth
done

#Average depth for aligned regions
for i in *.depth ; do
 	n=$(awk '{sum+=$3} END { print "Average = ",sum/NR}' $i)
 	echo ${i%.depth} $n >> avg.depth.alignedRegions.txt 
done

#make depth files for just exon regions, so I can look at average depth and not just number of positions aligning
for i in *.bam ; do 
	awk -v OFS='\t' '{print "Chr"$0}' ../ref_Sb313.cds.bed | samtools depth $i -b - > ${i%.bam}.depth.refexons
done

#Average depth for exon regions
for i in *.depth.refexons ; do 
	n=$(awk '{sum+=$3} END { print "Average = ",sum/NR}' $i)
 	echo ${i%.depth.refexons} $n >> avg.depth.refExons.txt
done

sed 's/Average = //g' avg.depth.refExons.txt | sed 's/_anchorwave//g' - | sed 's/.maf//g' - | sed 's/Sb313_//g' - | tr " " "\t" >avg.depth.refExons.tsv
sed 's/Average = //g' avg.depth.alignedRegions.txt | sed 's/_anchorwave//g' - | sed 's/.maf//g' - | sed 's/Sb313_//g' - | tr " " "\t" >avg.depth.alignedRegions.tsv
```

### Option 3: use MAF > PSL > BED to retrieve all aligned regions between ref and query

### Option 4: use WGAbed to convert MAF to Bed (Slow)

### Option 5: Look at the number of anchors created for each genome relative to the number of genes in the reference gff
```
awk -v OFS='\t' '$3 ~ /gene/ {print $0}' Sbicolor_313/Sbicolor_313_v3.1/Sbicolor_313_v3.1.gene.gff3 | wc -l 
#number of genes in ref gff3
#34211

for i in *.anchors ; do 
	a=$(grep -c "Sobic" $i)
	l=$(grep -c "localAlignment" $i)
	echo $i $a $l >> anchor.counts
done
sed -i 's/ /\t/g' anchor.counts
sed -i 's/Sb313_//g' anchor.counts
sed -i 's/_anchorwave.anchors//g' anchor.counts
```

## Make  the overall fractionation calls

2. Differentiate between aligned regions and unaligned regions
a. which reference exons are in the aligned regions? 
b. which reference exons are in the unaligned regions? --> assign NA

3. For the aligned regions with ref exons, which reference exons intersect with an indel?
a. in aligned region, no deletion intersecting --> retained
b. in aligned region, deletion intersecting --> fractionated
	i. Any way to count how many base pairs are in the deletion from the exon? To get at proportion? 
		Probably better to leave that for the gene level since deletions have fuzzy boundaries

_Make sure GVCFs are moved out of the directory used for running the VCF metrics plugin_
		
```edited from slurm_08.00.gvcfDelIntersect.sh

#split reference exons by reference CHR

for j in 01 02 03 04 05 06 07 08 09 10 ; do 
	awk -v OFS='\t' -v cn=$(echo $j) '$1 == cn {print $0}' /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed > /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${j}.bed
done

ml bedtools2

#you have to make the GVCF into a bed file with a start and end position (either with the END= info for non-variant entries or the length of the ref allele for variant entries). Otherwise the intersect function will only intersect the
 points and not the actual blocks within the GVCF
 
for i in GVCF/Sb313*.gvcf.gz; do
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	o=${i#GVCF/} #Get the GVCF name
	zcat $i | grep -v "^#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - | \ #make it into an actual bed file with start and end positions (otherwise the intersect will just be intersecting the points and not actual blocks)
	sed 's/END=//g' - | sed 's/ASM_Strand=//g' - | \
	bedtools intersect -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | \
	cut -f 4 | sort | uniq > ${o%.gvcf.gz}.aligned.IDs.txt
	zcat $i | grep -v "^#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - | \
	sed 's/END=//g' - | sed 's/ASM_Strand=//g' -| \
	bedtools intersect -v -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | \
	cut -f 4 | sort | uniq > ${o%.gvcf.gz}.noalignment.IDs.txt
done

##How to check if the above worked
for i in *align*.IDs.txt ; do 
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' - ) #get ref chr
	c=$(wc -l $i)
	o=$(echo $i | cut -f 1-2 -d ".") #stem
	echo $o $n $c ${i%.IDs.txt} >> QC.align.calls.txt #stem-name, ref Chr, count, file name (aligned or not)
done

sed -i 's/ /\t/g' QC.align.calls.txt

#numbers come from wc -l ref_Sb313.cds.${refchr}.bed

awk -v OFS='\t' '{if($2 == "01") print $0,13388 ; else print $0}' QC.align.calls.txt > QC.align.calls.txt.1
awk -v OFS='\t' '{if($2 == "02") print $0,8683 ; else print $0}' QC.align.calls.txt.1 > QC.align.calls.txt.2
awk -v OFS='\t' '{if($2 == "03") print $0,10580 ; else print $0}' QC.align.calls.txt.2 > QC.align.calls.txt.3
awk -v OFS='\t' '{if($2 == "04") print $0,9358 ; else print $0}' QC.align.calls.txt.3 > QC.align.calls.txt.4
awk -v OFS='\t' '{if($2 == "05") print $0,2311 ; else print $0}' QC.align.calls.txt.4 > QC.align.calls.txt.5
awk -v OFS='\t' '{if($2 == "06") print $0,6169 ; else print $0}' QC.align.calls.txt.5 > QC.align.calls.txt.6
awk -v OFS='\t' '{if($2 == "07") print $0,4331 ; else print $0}' QC.align.calls.txt.6 > QC.align.calls.txt.7
awk -v OFS='\t' '{if($2 == "08") print $0,2898 ; else print $0}' QC.align.calls.txt.7 > QC.align.calls.txt.8
awk -v OFS='\t' '{if($2 == "09") print $0,5743 ; else print $0}' QC.align.calls.txt.8 > QC.align.calls.txt.9
awk -v OFS='\t' '{if($2 == "10") print $0,5808 ; else print $0}' QC.align.calls.txt.9 > QC.align.calls.txt

rm QC.align.calls.txt.[1-9]

awk -v OFS='\t' '{if ($1 == last) {if(prevCnt+$3 == $6) print $0,"TRUE" ; else print $0, "FALSE"}; {last=$1; prevCnt=$3}}' QC.align.calls.txt | grep "FALSE"

#Got all TRUEs (700 as expected)

#finds deletions genotyped in the gvcfs
#Ignores if there is a deletion but only N as the ALT
#Makes the size of the deletion the ref allele - the alt allele lengths

for i in GVCF/Sb313*.gvcf.gz ; do 
o=${i#GVCF/} #Get the GVCF name
zcat $i | grep -v "#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{print $1,$2,$4,$5,$13}' |
	sed "s/<NON_REF>/N/g" | \
	awk -v OFS="\t" '{split($4,a,/,/); for(i = 1; i <= length(a); ++i) if(a[i] != "N" && length(a[i]) < length($3)) print $1,$2,$3,$4,a[i],length($3)-length(a[i]),$5}' |
	awk -v OFS="\t" '{print $1,$2,$2+$6,$3,$4,$5,$6,$7}' | \
	tr ":" "\t" | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${o%.gvcf.gz}.dels.bed
done

#fields are: SBChr, Sbstart, Sbstop, REF, ALT, SingularALT, REF_length (AKA deletion size), genotype

for i in *.dels.bed ; do
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	id=$(echo $i | sed 's/.dels.bed/.aligned.IDs.txt/g')
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -a - -b $i > ${i%.dels.bed}.refExons.deleted
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -v -a - -b $i > ${i%.dels.bed}.refExons.retained
	
	cut -f 4 ${i%.dels.bed}.refExons.deleted | sort | uniq > ${i%.dels.bed}.refExons.deleted.IDs.txt
	cut -f 4 ${i%.dels.bed}.refExons.retained | sort | uniq > ${i%.dels.bed}.refExons.retained.IDs.txt
done

####MAKE SURE THIS IS DONE IN A CLEAN DIRECTORY SO THERE'S NO OVERLAP
for i in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
 echo ${i}.bin1.noAlignment_CDS_ID > ${i}.bin1.allchr.refExons.noalign
 echo ${i}.bin2.noAlignment_CDS_ID > ${i}.bin2.allchr.refExons.noalign
 echo ${i}.bin1.deleted_CDS_ID > ${i}.bin1.allchr.refExons.deleted
 echo ${i}.bin2.deleted_CDS_ID > ${i}.bin2.allchr.refExons.deleted
 echo ${i}.bin1.retained_CDS_ID > ${i}.bin1.allchr.refExons.retained
 echo ${i}.bin2.retained_CDS_ID > ${i}.bin2.allchr.refExons.retained
 for j in 01 02 03 04 05 06 07 08 09 10; do 
 	cat Sb313_${i}_Chr${j}.bin1.noalignment.IDs.txt >> ${i}.bin1.allchr.refExons.noalign
 	cat Sb313_${i}_Chr${j}.bin2.noalignment.IDs.txt >> ${i}.bin2.allchr.refExons.noalign
 	cat Sb313_${i}_Chr${j}.bin1.refExons.deleted.IDs.txt >> ${i}.bin1.allchr.refExons.deleted
 	cat Sb313_${i}_Chr${j}.bin2.refExons.deleted.IDs.txt >> ${i}.bin2.allchr.refExons.deleted
 	cat Sb313_${i}_Chr${j}.bin1.refExons.retained.IDs.txt  >> ${i}.bin1.allchr.refExons.retained
 	cat Sb313_${i}_Chr${j}.bin2.refExons.retained.IDs.txt  >> ${i}.bin2.allchr.refExons.retained
 done
done	

# QC check
for i in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	for b in bin1 bin2 ; do
		a=$(wc -l ${i}.${b}.allchr.refExons.noalign)
		d=$(wc -l ${i}.${b}.allchr.refExons.deleted)
		r=$(wc -l ${i}.${b}.allchr.refExons.retained)
		echo $i $b $a $d $r >> QC.totalcalls.txt
	done
done

sed -i 's/ /\t/g' QC.totalcalls.txt 

#Recall that there were some ref exons that were on scaffolds, not chr, so only looking at those on chr
#numbers come from wc -l ../../ref_Sb313.cds.[0-9]*.bed 
#look at total then add 3 for the headers
69269 + 3

awk -v OFS="\t" '{if ($3+$5+$7 == 69272) print $0, "TRUE" ; else print $0, "FALSE"}' QC.totalcalls.txt > QC.totalcalls.tsv
#Got all trues and 70 of them like expected

for i in TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	for b in bin1 bin2 ; do
		for j in 01 02 03 04 05 06 07 08 09 10 ; do
		a=$(wc -l Sb313_${i}_Chr${j}.${b}.noalignment.IDs.txt)
		d=$(wc -l Sb313_${i}_Chr${j}.${b}.refExons.deleted.IDs.txt)
		r=$(wc -l Sb313_${i}_Chr${j}.${b}.refExons.retained.IDs.txt)
		echo ${i} ${b} ${j} ${a} ${d} ${r} >> QC.totalcalls.byRefChr.txt
		done
	done
done
sed 's/ /\t/g' QC.totalcalls.byRefChr.txt > QC.totalcalls.byRefChr.tsv
```

*From Maggie's manual curation (see QC scripts), large deletions and large/segregating deletions tended to be false AW calls*
Turned out to be due to the strandedness in changing deletions to bed file intervals
Edited 2/23 and rerun

For the TdKS and 4to1 Zn Zd run (so I don't have to redo everything):
```
ml bedtools2

#you have to make the GVCF into a bed file with a start and end position (either with the END= info for non-variant entries or the length of the ref allele for variant entries). Otherwise the intersect function will only intersect the
 points and not the actual blocks within the GVCF
 
for i in GVCF/Sb313*4to1*.gvcf.gz ; do
	n=$(echo $i | cut -f 4 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	o=${i#GVCF/} #Get the GVCF name
	zcat $i | grep -v "^#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - | \ #make it into an actual bed file with start and end positions (otherwise the intersect will just be intersecting the points and not actual blocks)
	sed 's/END=//g' - | sed 's/ASM_Strand=//g' - | \
	bedtools intersect -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | \
	cut -f 4 | sort | uniq > ${o%.gvcf.gz}.aligned.IDs.txt
	zcat $i | grep -v "^#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - | \
	sed 's/END=//g' - | sed 's/ASM_Strand=//g' -| \
	bedtools intersect -v -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | \
	cut -f 4 | sort | uniq > ${o%.gvcf.gz}.noalignment.IDs.txt
done
for i in GVCF/Sb313*TdKS*.gvcf.gz ; do
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	o=${i#GVCF/} #Get the GVCF name
	zcat $i | grep -v "^#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - | \ #make it into an actual bed file with start and end positions (otherwise the intersect will just be intersecting the points and not actual blocks)
	sed 's/END=//g' - | sed 's/ASM_Strand=//g' - | \
	bedtools intersect -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | \
	cut -f 4 | sort | uniq > ${o%.gvcf.gz}.aligned.IDs.txt
	zcat $i | grep -v "^#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{if($12 ~/END/) print $1,$2,$12,$4,$5,$11,$14 ; else print $1,$2,$2+length($4),$4,$5,$11,$13}' - | \
	sed 's/END=//g' - | sed 's/ASM_Strand=//g' -| \
	bedtools intersect -v -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed -b - | \
	cut -f 4 | sort | uniq > ${o%.gvcf.gz}.noalignment.IDs.txt
done


##How to check if the above worked
for i in *align*.IDs.txt ; do 
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' - ) #get ref chr
	c=$(wc -l $i)
	o=$(echo $i | cut -f 1-2 -d ".") #stem
	echo $o $n $c ${i%.IDs.txt} >> QC.align.calls.txt #stem-name, ref Chr, count, file name (aligned or not)
done

sed -i 's/ /\t/g' QC.align.calls.txt

#numbers come from wc -l ref_Sb313.cds.${refchr}.bed

awk -v OFS='\t' '{if($2 == "01") print $0,13388 ; else print $0}' QC.align.calls.txt > QC.align.calls.txt.1
awk -v OFS='\t' '{if($2 == "02") print $0,8683 ; else print $0}' QC.align.calls.txt.1 > QC.align.calls.txt.2
awk -v OFS='\t' '{if($2 == "03") print $0,10580 ; else print $0}' QC.align.calls.txt.2 > QC.align.calls.txt.3
awk -v OFS='\t' '{if($2 == "04") print $0,9358 ; else print $0}' QC.align.calls.txt.3 > QC.align.calls.txt.4
awk -v OFS='\t' '{if($2 == "05") print $0,2311 ; else print $0}' QC.align.calls.txt.4 > QC.align.calls.txt.5
awk -v OFS='\t' '{if($2 == "06") print $0,6169 ; else print $0}' QC.align.calls.txt.5 > QC.align.calls.txt.6
awk -v OFS='\t' '{if($2 == "07") print $0,4331 ; else print $0}' QC.align.calls.txt.6 > QC.align.calls.txt.7
awk -v OFS='\t' '{if($2 == "08") print $0,2898 ; else print $0}' QC.align.calls.txt.7 > QC.align.calls.txt.8
awk -v OFS='\t' '{if($2 == "09") print $0,5743 ; else print $0}' QC.align.calls.txt.8 > QC.align.calls.txt.9
awk -v OFS='\t' '{if($2 == "10") print $0,5808 ; else print $0}' QC.align.calls.txt.9 > QC.align.calls.txt

rm QC.align.calls.txt.[1-9]

awk -v OFS='\t' '{if ($1 == last) {if(prevCnt+$3 == $6) print $0,"TRUE" ; else print $0, "FALSE"}; {last=$1; prevCnt=$3}}' QC.align.calls.txt | grep "FALSE"

#Got all TRUEs (700 as expected)

#finds deletions genotyped in the gvcfs
#Ignores if there is a deletion but only N as the ALT
#Makes the size of the deletion the ref allele - the alt allele lengths

for i in GVCF/Sb313*4to1*.gvcf.gz GVCF/Sb313*TdKS*.gvcf.gz ; do 
o=${i#GVCF/} #Get the GVCF name
zcat $i | grep -v "#" - | tr ';' '\t' | \
	awk -v OFS="\t" '{print $1,$2,$4,$5,$13}' |
	sed "s/<NON_REF>/N/g" | \
	awk -v OFS="\t" '{split($4,a,/,/); for(i = 1; i <= length(a); ++i) if(a[i] != "N" && length(a[i]) < length($3)) print $1,$2,$3,$4,a[i],length($3)-length(a[i]),$5}' |
	awk -v OFS="\t" '{print $1,$2,$2+$6,$3,$4,$5,$6,$7}' | \
	tr ":" "\t" | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8}' > ${o%.gvcf.gz}.dels.bed
done


#fields are: SBChr, Sbstart, Sbstop, REF, ALT, SingularALT, REF_length (AKA deletion size), genotype

for i in *TdKS*.dels.bed ; do
	n=$(echo $i | cut -f 3 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	id=$(echo $i | sed 's/.dels.bed/.aligned.IDs.txt/g')
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -a - -b $i > ${i%.dels.bed}.refExons.deleted
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -v -a - -b $i > ${i%.dels.bed}.refExons.retained
	cut -f 4 ${i%.dels.bed}.refExons.deleted | sort | uniq > ${i%.dels.bed}.refExons.deleted.IDs.txt
	cut -f 4 ${i%.dels.bed}.refExons.retained | sort | uniq > ${i%.dels.bed}.refExons.retained.IDs.txt
done
for i in *4to1*.dels.bed ; do
	n=$(echo $i | cut -f 4 -d "_" | cut -f 1 -d "." | sed 's/Chr//g' -) #get the ref chr
	id=$(echo $i | sed 's/.dels.bed/.aligned.IDs.txt/g')
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -a - -b $i > ${i%.dels.bed}.refExons.deleted
	grep -f $id /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${n}.bed | bedtools intersect -v -a - -b $i > ${i%.dels.bed}.refExons.retained
	
	cut -f 4 ${i%.dels.bed}.refExons.deleted | sort | uniq > ${i%.dels.bed}.refExons.deleted.IDs.txt
	cut -f 4 ${i%.dels.bed}.refExons.retained | sort | uniq > ${i%.dels.bed}.refExons.retained.IDs.txt
done

####MAKE SURE THIS IS DONE IN A CLEAN DIRECTORY SO THERE'S NO OVERLAP
for i in TdKS ZdGigi_4to1 ZdMomo_4to1 ZnPI615697_4to1 ; do 
 echo ${i}.bin1.noAlignment_CDS_ID > ${i}.bin1.allchr.refExons.noalign
 echo ${i}.bin2.noAlignment_CDS_ID > ${i}.bin2.allchr.refExons.noalign
 echo ${i}.bin1.deleted_CDS_ID > ${i}.bin1.allchr.refExons.deleted
 echo ${i}.bin2.deleted_CDS_ID > ${i}.bin2.allchr.refExons.deleted
 echo ${i}.bin1.retained_CDS_ID > ${i}.bin1.allchr.refExons.retained
 echo ${i}.bin2.retained_CDS_ID > ${i}.bin2.allchr.refExons.retained
 for j in 01 02 03 04 05 06 07 08 09 10; do 
 	cat Sb313_${i}_Chr${j}.bin1.noalignment.IDs.txt >> ${i}.bin1.allchr.refExons.noalign
 	cat Sb313_${i}_Chr${j}.bin2.noalignment.IDs.txt >> ${i}.bin2.allchr.refExons.noalign
 	cat Sb313_${i}_Chr${j}.bin1.refExons.deleted.IDs.txt >> ${i}.bin1.allchr.refExons.deleted
 	cat Sb313_${i}_Chr${j}.bin2.refExons.deleted.IDs.txt >> ${i}.bin2.allchr.refExons.deleted
 	cat Sb313_${i}_Chr${j}.bin1.refExons.retained.IDs.txt  >> ${i}.bin1.allchr.refExons.retained
 	cat Sb313_${i}_Chr${j}.bin2.refExons.retained.IDs.txt  >> ${i}.bin2.allchr.refExons.retained
 done
done	

# QC check
for i in TdKS ZdGigi_4to1 ZdMomo_4to1 ZnPI615697_4to1  TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	for b in bin1 bin2 ; do
		a=$(wc -l ${i}.${b}.allchr.refExons.noalign)
		d=$(wc -l ${i}.${b}.allchr.refExons.deleted)
		r=$(wc -l ${i}.${b}.allchr.refExons.retained)
		echo $i $b $a $d $r >> QC.totalcalls.txt
	done
done

sed -i 's/ /\t/g' QC.totalcalls.txt 

#Recall that there were some ref exons that were on scaffolds, not chr, so only looking at those on chr
#numbers come from wc -l ../../ref_Sb313.cds.[0-9]*.bed 
#look at total then add 3 for the headers
69269 + 3

awk -v OFS="\t" '{if ($3+$5+$7 == 69272) print $0, "TRUE" ; else print $0, "FALSE"}' QC.totalcalls.txt > QC.totalcalls.tsv
#Got all trues and 70 of them like expected

for i in TdKS ZdGigi_4to1 ZdMomo_4to1 ZnPI615697_4to1  TdFL ZdGigi ZdMomo ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZnPI615697 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do 
	for b in bin1 bin2 ; do
		for j in 01 02 03 04 05 06 07 08 09 10 ; do
		a=$(wc -l Sb313_${i}_Chr${j}.${b}.noalignment.IDs.txt)
		d=$(wc -l Sb313_${i}_Chr${j}.${b}.refExons.deleted.IDs.txt)
		r=$(wc -l Sb313_${i}_Chr${j}.${b}.refExons.retained.IDs.txt)
		echo ${i} ${b} ${j} ${a} ${d} ${r} >> QC.totalcalls.byRefChr.txt
		done
	done
done
sed 's/ /\t/g' QC.totalcalls.byRefChr.txt > QC.totalcalls.byRefChr.tsv
```


4. Redo the analysis with the new calls

```
1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8787:nova18-26:34033 snodgras@nova.its.iastate.edu

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: snodgras
   password: yWKb1zJrR7XMcdCttavM

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      scancel -f 5711897

```


## To Use the VCFs for identifying multiple origins from the GVCFs
Need to trim out the too big deletions
Then combine GVCFs into VCFS
Then filter to only indels
Convert VCFs to bed-like files

Compare:
1. exonic vs. non-exonic
2. single deletions/biallelic vs. nested
3. within exon boundaries vs. extending past


## Trimming the GVCFs
(use the `test_mafs` first)
```
for i in *gvcf.gz ; do 
echo module load samtools \; module load gatk \; gatk LeftAlignAndTrimVariants -O trimmed_${i} -V ${i} -R /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.fasta --max-indel-length 9101263 &> gf_stdout.txt >> splitmaf.trimGVCF.cmds.txt
done
cp ../../../scripts/makeSLURM.py .
#edit makeSLURM.py to be only 1 hour of wall time
python makeSLURM.py 1 splitmaf.trimGVCF.cmds.txt

for i in *gvcf.gz ; do echo $i >> job.key.tmp ; done
awk -v OFS='\t' '{print $0,NR-1}' job.key.tmp > job.key
rm job.key.tmp

while read -r id job ; do 
	sed -i "s/splitmaf/${id}/g" splitmaf.trimGVCF.cmds_"$job".sub
done < job.key

sbatch splitmaf.trimGVCF.cmds_0.sub
for i in {1..7} ; do sbatch --dependency=afterok:4858025 splitmaf.trimGVCF.cmds_${i}.sub ;done
```

Running on the full set of genomes
```
for i in *gvcf.gz ; do 
echo module load samtools \; module load gatk \; gatk LeftAlignAndTrimVariants -O trimmed_${i} -V ${i} -R /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.fasta --max-indel-length 9101263 &> gf_stdout.txt >> splitmaf.trimGVCF.cmds.txt
done
python ../../scripts/makeSLURM.py 1 splitmaf.trimGVCF.cmds.txt

for i in *gvcf.gz ; do echo $i >> job.key.tmp ; done
awk -v OFS='\t' '{print $0,NR-1}' job.key.tmp > job.key
rm job.key.tmp

while read -r id job ; do 
	sed -i "s/splitmaf/${id}/g" splitmaf.trimGVCF.cmds_"$job".sub
done < job.key

for i in splitmaf.trimGVCF.cmds*sub ; do sbatch $i ; done
```
TdKS and ZnZd 4to1 GVCFs
_MUST RUN IndexFeatureFile on input to get the trimmed gvcf_
```
for i in GVCF/*TdKS* GVCF/*4to1* ; do
echo module load samtools \; module load gatk \; gatk IndexFeatureFile -I ${i} \; gatk LeftAlignAndTrimVariants -O trimmed_${i#GVCF/} -V ${i} -R /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.fasta --max-indel-length 9101263 &> gf_stdout.txt >> splitmaf.TdKSZnZd.trimGVCF.cmds.txt
done
python ../../scripts/makeSLURM.py 1 splitmaf.TdKSZnZd.trimGVCF.cmds.txt

for i in GVCF/*TdKS* GVCF/*4to1* ; do echo $i >> TdKSZnZd.job.key.tmp ; done
awk -v OFS='\t' '{print $0,NR-1}' TdKSZnZd.job.key.tmp > TdKSZnZd.job.key
rm TdKSZnZd.job.key.tmp

while read -r id job ; do 
	sed -i "s/splitmaf/${id#GVCF/}/g" splitmaf.TdKSZnZd.trimGVCF.cmds_"$job".sub
done < TdKSZnZd.job.key

for i in splitmaf.TdKSZnZd.trimGVCF.cmds_*.sub ; do sbatch $i ; done
```
Didn't work for `Zn_4to1`
So for it: 

```
ml bcftools gatk samtools
 bcftools sort GVCF/Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz -o GVCF/sorted_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz
 cd GVCF/
 mv sorted_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz sorted_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf
 bgzip sorted_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf 
 cd ..
 gatk IndexFeatureFile -I GVCF/sorted_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz 
 gatk CreateSequenceDictionary REFERENCE=/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.fasta OUTPUT=/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.dict
 gatk LeftAlignAndTrimVariants -O trimmed_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz -V GVCF/sorted_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz  -R /work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.clean.fasta --max-indel-length 9101263 &> gf_stdout.txt
awk '$0 ~ /Indel is too long/  { print }' gf_stdout.txt |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
vcftools --gzvcf trimmed_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz --exclude-positions positions.txt --recode --recode-INFO-all --out ready_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz.recode.vcf
 gatk IndexFeatureFile -I ready_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz.recode.vcf
```

## Then manually remove the too large indels

(use the `test_mafs` first)
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

for i in Sb313*.gvcf.gz ; do 
	awk '$0 ~ /Indel is too long/  { print }' ${i}.trimGVCF.cmds_*.e* |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
	vcftools --gzvcf trimmed_${i} --exclude-positions positions.txt --recode --recode-INFO-all --out ready_${i}
done

ml gatk
for i in ready*recode.vcf ; do gatk IndexFeatureFile -I $i ;done


for i in GVCF/*TdKS*.gvcf.gz GVCF/*4to1*.gvcf.gz ; do 
	awk '$0 ~ /Indel is too long/  { print }' ${i#GVCF/}.TdKSZnZd.trimGVCF.cmds_*.e* |cut -f 1 -d ";" |awk '{print $NF}' |sed 's/:/\t/g' > positions.txt
	vcftools --gzvcf trimmed_${i#GVCF/} --exclude-positions positions.txt --recode --recode-INFO-all --out ready_${i#GVCF/}
done

for i in ready*TdKS*recode.vcf ready*4to1*recode.vcf ; do gatk IndexFeatureFile -I $i ; done
```
output is like: `ready_Sb313_TdFL_Chr10.bin1.gvcf.gz.recode.vcf`

Some of the large indels remain
```slurm-04.6.removeremaininglargeindels.sh
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=8:00:00

for i in ready_Sb313_*Chr01*bin1*vcf ; do
	echo $i >> problem_indels.txt
	awk -v OFS='\t' 'length($4) > 9101263 || length($5) - 10 > 9101263 {print $2, "ref length:"length($4), "alt length:"length($5)-10}' ${i} >> problem_indels.txt
done

for j in 02 03 04 05 06 07 08 09 10 ; 
	do for i in ready_Sb313_*Chr${j}*bin1*.vcf ; 
		do echo $i >> problem_indels.txt ; 
		awk -v OFS='\t' 'length($4) > 9101263 || length($5) - 10 > 9101263 {print $2, "ref length:"length($4), "alt length:"length($5)-10}' ${i} >> problem_indels.txt ; 
	done ; 
done  

for j in 01 02 03 04 05 06 07 08 09 10 ; 
	do for i in ready_Sb313_*Chr${j}*bin2*.vcf ; 
		do echo $i >> problem_indels.txt ; 
	    awk -v OFS='\t' 'length($4) > 9101263 || length($5) - 10 > 9101263 {print $2, "ref length:"length($4), "alt length:"length($5)-10}' ${i} >> problem_indels.txt ; 
	done ; 
done  

#only Chr01.bin1 files had this issue
#ZmCML333 ZmCML52 ZmIL14H ZmTx303 ZvTIL11
for i in ready_Sb313_*_Chr01.bin1.gvcf.gz.recode.vcf ; do 
	awk -v OFS='\t' 'length($4) < 9101263 && length($5) - 10 < 9101263 {print $0}' ${i} > pass_${i}
done
sed -i 's/ready/pass_ready/g' /work/LAS/mhufford-lab/snodgras/Fractionation/scripts/slurm_06.1.splitmaf.bin1.gvcf2vcf.sh
```

## Combine the GVCFs into a VCF
`06.{chr}.splitmaf.{bin}.gvcf2vcf.sh`
```
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=48:00:00
#SBATCH --job-name=gatk-chr1
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

ml gatk
ref=/work/LAS/mhufford-lab/snodgras/Fractionation/Sb313.chr1.fasta
ml samtools
# index
samtools faidx $ref
# dict
gatk CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${ref%.*}.dict
# db import
gatk --java-options "-Xmx128g -Xms5g" GenomicsDBImport \
-V ready_Sb313_TdFL_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_TdKS_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdGigi_4to1_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZdMomo_4to1_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZhRIMHU001_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB73_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmB97_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML103_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML228_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML247_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML277_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML322_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML333_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML52_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmCML69_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmHP301_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmIL14H_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmKi11_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmKi3_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmKy21_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmM162W_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmM37W_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmMo18W_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmMS71_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmNC350_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmNC358_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmOh43_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmOh7b_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmP39_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmTx303_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZmTzi8_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZnPI615697_4to1_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL01_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZvTIL11_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZxTIL18_Chr01.bin1.gvcf.gz.recode.vcf \
-V ready_Sb313_ZxTIL25_Chr01.bin1.gvcf.gz.recode.vcf \
--batch-size 1 \
--genomicsdb-workspace-path chr1.bin1_gatkDBimport \
-L 01 \
--genomicsdb-segment-size 1048576000 --genomicsdb-vcf-buffer-size 10000000000 
#--tmp-dir $TMPDIR

# vcf output
gatk --java-options "-Xmx50g" GenotypeGVCFs \
-R $ref \
-V gendb://chr1.bin1_gatkDBimport \
--cloud-prefetch-buffer 10000 --cloud-index-prefetch-buffer 10000 --genomicsdb-max-alternate-alleles 110 --max-alternate-alleles 100 --gcs-max-retries 1000 \
-O chr1.bin1.vcf
```
To make copies of the above script for each chromosome and bin:

```
cp slurm_06.1.splitmaf.bin1.gvcf2vcf.sh slurm_06.1.splitmaf.bin2.gvcf2vcf.sh 
sed -i 's/bin1/bin2/g' slurm_06.1.splitmaf.bin2.gvcf2vcf.sh 
#manually comment out the dictionary making lines of the bin2 script

for i in {2..10} ; do 
	for b in bin1 bin2 ; do 
		cp slurm_06.1.splitmaf.${b}.gvcf2vcf.sh slurm_06.${i}.splitmaf.${b}.gvcf2vcf.sh 
	done
done

for c in slurm_06.2.splitmaf*sh ; do sed -i 's/01/02/g' $c ; sed -i 's/chr1/chr2/g' $c ; done
for c in slurm_06.3.splitmaf*sh ; do sed -i 's/01/03/g' $c ; sed -i 's/chr1/chr3/g' $c ; done
for c in slurm_06.4.splitmaf*sh ; do sed -i 's/01/04/g' $c ; sed -i 's/chr1/chr4/g' $c ; done
for c in slurm_06.5.splitmaf*sh ; do sed -i 's/01/05/g' $c ; sed -i 's/chr1/chr5/g' $c ; done
for c in slurm_06.6.splitmaf*sh ; do sed -i 's/01/06/g' $c ; sed -i 's/chr1/chr6/g' $c ; done
for c in slurm_06.7.splitmaf*sh ; do sed -i 's/01/07/g' $c ; sed -i 's/chr1/chr7/g' $c ; done
for c in slurm_06.8.splitmaf*sh ; do sed -i 's/01/08/g' $c ; sed -i 's/chr1/chr8/g' $c ; done
for c in slurm_06.9.splitmaf*sh ; do sed -i 's/01/09/g' $c ; sed -i 's/chr1/chr9/g' $c ; done
for c in slurm_06.10.splitmaf*sh ; do sed -i 's/01/10/g' $c ; sed -i 's/chr1/chr10/g' $c ; done

#go back and change ZhRIMHU00[2-9] to ZhRIMHU001
#go back and change HP30[2-9] to HP301
#go back and change ZvTIL0[2-9] to ZvTIL01
#make sure there are no `chr1.bin1_gatkDBimport` pre-existing before running
```
Chr1 bin1 continues to run into buffer overflow issues (even after extra pruning and trimming)
Chr10 bin1 has an issue with the 32nd file (ready_Sb313_ZnPI615697_4to1_Chr10.bin1.gvcf.gz.recode.vcf):
 what():  VCF2TileDBException : Incorrect cell order found - cells must be in column major order. Previous cell: [ 31, 
681363017 ] current cell: [ 31, 681363017 ].
The most likely cause is unexpected data in the input file:
(a) A VCF file has two lines with the same genomic position
(b) An unsorted CSV file
(c) Malformed VCF file (or malformed index)

_Since Zn won't be used in analyses, I'm trying to run the script again just commenting out the line with Zn_
If it works, we'll just use that rather than trying to troubleshoot the above. 

Still failed
troubleshooting in `troubleshoot_gvcf2vcf_error/` directory
Trying: 1. get rid of Zn as an option, 2. subset to only 10 genomes, 3. exclude the region that throws the error
Need to index the gvcfs before trying to run the above options by doing:
`ml gatk ; for i in pass*vcf ; do gatk IndexFeatureFile -I $i ;done`

Sorghum Chr1 is 80874895 bp long for reference

the exclude the 40M run didn't work
the run without Zn worked

Filtering out Zn from all the other vcfs and putting in a separate directory for further work
```
mkdir noZn_vcf
mv *noZn*.vcf noZn_vcf/.

ml bcftools

for c in {2..9} ; do 
	for b in 1 2 ; do 
		s=$(echo Sb313_ZnPI615697_4to1_Chr0${c}.bin${b})
		bcftools view -s ^${s} -o noZn_vcf/chr0${c}.bin${b}.noZn.vcf chr${c}.bin${b}.vcf
	done
done

```


## Filter VCFs to indels only
`slurm_07.{chr}.splitmaf.subsetvcf.sh`
```
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

gatk SelectVariants -V chr1.bin1.vcf -O chr1.bin1.indelonly.vcf --select-type-to-include INDEL
gatk SelectVariants -V chr1.bin2.vcf -O chr1.bin2.indelonly.vcf --select-type-to-include INDEL
```

## Reformat VCFs to bed-like 
```
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=1:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#ml gatk
#gatk VariantsToTable -V chr10.indelonly.vcf -O chr10.indelonly.table -F CHROM -F POS -F TYPE -F NCALLED -GF GT -GF GQ

ml bcftools
bcftools annotate -x INFO,^FORMAT/GT chr10.bin1.indelonly.vcf  > chr10.bin1.indelonly_reformatted.vcf
bcftools annotate -x INFO,^FORMAT/GT chr10.bin2.indelonly.vcf  > chr10.bin2.indelonly_reformatted.vcf
```
Lump both together to get
```
ml bcftools
for c in 01 02 03 04 05 06 07 08 09 10 ; do 
	for b in 1 2 ; do
		bcftools view --types indels chr${c}.bin${b}.noZn.vcf | bcftools annotate -x INFO,^FORMAT/GT - > chr${c}.bin${b}.indelonly_reformatted.vcf
	done
done
```

Separate out the deletion variants from insertion variants and convert to bed

```
for i in *.indelonly_reformatted.vcf ; do 
	if [ ! -f ${i%.vcf}_headerless.vcf ] ; then
		grep -v "^##" ${i} > ${i%.vcf}_headerless.vcf ;
	fi
	head -n 1 ${i%.vcf}_headerless.vcf > ${i%.indelonly_reformatted_headerless.vcf}.delonly.vcf
	cat ${i%.vcf}_headerless.vcf | awk -v OFS='\t' '{split($5,a,/,/); if(length(a[1]) < length($4) || (length(a[2]) < length($4) && a[2] != "")) print $0}' - >> ${i%.indelonly_reformatted_headerless.vcf}.delonly.vcf
done

echo Done making deletion only vcfs!

for i in *delonly.vcf ; do 
	awk -v OFS='\t' '{split($5,a,/,/); 
		print $1,$2,$2+length($4),$4,$5,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45}' ${i} > ${i%.vcf}.bed
done

echo Done converting del only vcfs to bed files!
```
The names didn't come out right... :/
```
for c in 01 02 03 04 05 06 07 08 09 10 ; do for b in 1 2 ; do 
	mv chr${c}.bin${b}.indelonly_reformatted.vcf.delonly.vcf chr${c}.bin${b}.delonly.vcf
	mv chr${c}.bin${b}.indelonly_reformatted.vcf.delonly.bed chr${c}.bin${b}.delonly.bed
	done ; done
```

Collapse variants (with the same start site)
```
# First check to see if there are any that share the same start position:
for i in *delonly.bed ; do echo $i ; sort -k 1,1 -k 2,2n $i | awk -v OFS='\t' '$2 == a { print $0 } {a=$2}' - | wc -l ; done
#none of them show multiple deletions starting at the same position
```

Collapse variants that overlap reciprocally by 80%
```
# remove the tab at the end of the line 
for i in *delonly.bed ; do sed -i 's/\t$//g' $i ; done

#subset deletions to just be those overlapping exons
for c in 01 02 03 04 05 06 07 08 09 10 ; do for b in bin1 bin2 ; do 
	bedtools intersect -wa  -a chr${c}.${b}.delonly.bed -b /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${c}.bed | uniq > chr${c}.${b}.del.exonic.bed
	done ; done

#keep in mind that this doesn't write which exons the deletions overlap; this is just to cut out the deletions that have NO overlap with exons
#keep in mind that when I did wc -l on all the *.exonic.bed, Chr5 and Chr10 had far fewer than the others

# create the reciprocal overlap between variants
for c in 01 02 03 04 05 06 07 08 09 10 ; do for b in bin1 bin2 ; do 
	bedtools intersect -a chr${c}.${b}.del.exonic.bed -b chr${c}.${b}.del.exonic.bed -f 0.8 -r -wa -wb > chr${c}.${b}.exonic.reciprocal80.bed
done ; done

#because it's reciprocal, it'll match the same variant to itself
#field 43 is the start position for the "overlap" and $44 is the stop position for the "overlap"

awk -v OFS='\t' 'if($2 != $43 && $3!=$44){print $1,$2,$3,$1":"$2"-"$3,$5,$42,$43,$44,$42":"$43"-"$44,$46}' chr10.bin1.exonic.reciprocal80.bed | sort -k2,2n -k7,7n | uniq > chr10.bin1.exonic.reciprocal80.dupremoved.key
bedtools groupby -i chr10.bin1.exonic.reciprocal80.dupremoved.key -grp 1-4 -c 9 -o collapse | awk -v OFS='\t' '$5 !~ p {print $0};{p=$1":"$2"-"$3}' - > chr10.bin1.exonic.collapse.key

# use awk to write the genotypes? 

awk -v OFS='\t' '{print $1":"$2"-"$3,$0}' chr10.bin1.del.exonic.bed > chr10.bin1.del.exonic.wIDs.bed

while read -r line ; do 
	id=$(echo $line | cut -f 1 -d " ") #create the id of the deletion
	if ! grep -q $id chr10.bin1.exonic.collapse.key #if the id IS NOT in the key of deletions that are collapsed
	then 
		echo $line | cut -f 2- -d " " >> chr10.bin1.del.exonic.collapsed.bed
	else
	if ! grep -q $id  chr10.bin1.del.exonic.collapsed.bed #if the id IS NOT already in the collapsed bed file (i.e. hasn't been added yet)  
		then
		grep $id chr10.bin1.exonic.collapse.key | cut -f 4-5 | tr "," "\t" | tr "\t" "\n" > grep.ids
		alt=$(grep $id chr10.bin1.exonic.collapse.key | cut -f 4-5 | tr "\t" ",")
		grep -w -f grep.ids chr10.bin1.del.exonic.wIDs.bed |
		awk -v OFS='\t' -v alt=$(echo $alt) '$3 < min {min=$3}; min == 0 {min=$3};
		$4 > max {max=$4}; max == 0 {max=$4};
		$8 != "." {g8 += $8}; 
		$9 != "." {g9 += $9} ; $10 != "." {g10 += $10} ;
		$11 != "." {g11 += $11} ; $12 != "." {g12 += $12} ; $13 != "." {g13 += $13} ;
		$14 != "." {g14 += $14} ; $15 != "." {g15 += $15} ; $16 != "." {g16 += $16} ;
		$17 != "." {g17 += $17} ; $18 != "." {g18 += $18} ; $19 != "." {g19 += $19} ;
		$20 != "." {g20 += $20} ; $21 != "." {g21 += $21} ; $22 != "." {g22 += $22} ;
		$23 != "." {g23 += $23} ; $24 != "." {g24 += $24} ; $25 != "." {g25 += $25} ;
		$26 != "." {g26 += $26} ; $27 != "." {g27 += $27} ; $28 != "." {g28 += $28} ;
		$29 != "." {g29 += $29} ; $30 != "." {g30 += $30} ; $31 != "." {g31 += $31} ;
		$32 != "." {g32 += $32} ; $33 != "." {g33 += $33} ; $14 != "." {g34 += $34} ;
		$35 != "." {g35 += $35} ; $36 != "." {g36 += $36} ; $37 != "." {g37 += $37} ;
		$38 != "." {g38 += $38} ; $39 != "." {g39 += $39} ; $40 != "." {g40 += $40} ;
		$41 != "." {g41 += $41}; $42 != "." {g42 += $42}
		END {print "10", min,max,"collapsedDEL",alt,"QUAL",g8,g9,g10,
		g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,
		g21,g22,g23,g24,g25,g26,g27,g28,g29,g30,
		g31,g32,g33,g34,g35,g36,g37,g38,g39,g40,
		g41,g42}' - | sed -E -e :1 -e 's/(^|\t)(\t|$)/\1\.\2/;t1' >> chr10.bin1.del.exonic.collapsed.bed
	fi
	fi
done < chr10.bin1.del.exonic.wIDs.bed

#testing this on the first 1000 records
while read -r line ; do id=$(echo $line | cut -f 1 -d " " ) ; if ! grep -q $id chr10.bin1.exonic.collapse.key  ; then echo $line | cut -f 2- -d " " >> chr10.bin1.del.exonic.collapsed.bed ; else if ! grep -q $id  chr10.bin1.del.exonic.collapsed.bed ; then grep $id chr10.bin1.exonic.collapse.key | cut -f 4-5 | tr "," "\t" | tr "\t" "\n" > grep.ids ; alt=$(grep $id chr10.bin1.exonic.collapse.key | cut -f 4-5 | tr "\t" ","); grep -w -f grep.ids chr10.bin1.del.exonic.wIDs.bed | awk -v OFS='\t' -v alt=$(echo $alt) '$3 < min {min=$3}; min == 0 {min=$3}; $4 > max {max=$4}; max == 0 {max=$4}; $8 != "." {g8 += $8}; $9 != "." {g9 += $9} ; $10 != "." {g10 += $10} ; $11 != "." {g11 += $11} ; $12 != "." {g12 += $12} ; $13 != "." {g13 += $13} ; $14 != "." {g14 += $14} ; $15 != "." {g15 += $15} ; $16 != "." {g16 += $16} ; $17 != "." {g17 += $17} ; $18 != "." {g18 += $18} ; $19 != "." {g19 += $19} ; $20 != "." {g20 += $20} ; $21 != "." {g21 += $21} ; $22 != "." {g22 += $22} ; $23 != "." {g23 += $23} ; $24 != "." {g24 += $24} ; $25 != "." {g25 += $25} ; $26 != "." {g26 += $26} ; $27 != "." {g27 += $27} ; $28 != "." {g28 += $28} ; $29 != "." {g29 += $29} ; $30 != "." {g30 += $30} ; $31 != "." {g31 += $31} ; $32 != "." {g32 += $32} ; $33 != "." {g33 += $33} ; $14 != "." {g34 += $34} ; $35 != "." {g35 += $35} ; $36 != "." {g36 += $36} ; $37 != "." {g37 += $37} ; $38 != "." {g38 += $38} ; $39 != "." {g39 += $39} ; $40 != "." {g40 += $40} ; $41 != "." {g41 += $41}; $42 != "." {g42 += $42} END {print "10", min,max,"collapsedDEL",alt,"QUAL",g8,g9,g10, g11,g12,g13,g14,g15,g16,g17,g18,g19,g20, g21,g22,g23,g24,g25,g26,g27,g28,g29,g30, g31,g32,g33,g34,g35,g36,g37,g38,g39,g40, g41,g42}' - |sed -E -e :1 -e 's/(^|\t)(\t|$)/\1\.\2/;t1' >> chr10.bin1.del.exonic.collapsed.bed ; fi ; fi ; done <test.bed

#not certain why I'm getting an awk cmd error but trying with the full file
# I think it works???

sed -i 's/ /\t/g' chr10.bin1.del.exonic.collapsed.bed 
```

Make the above code block a stand alone script
```collapseDels.sh
#!/bin/bash
inputbed=$1 #delonly.bed
chr=$2	#10
bin=$3	#bin1

ml bedtools2

#remove extra tabs from end of line
sed -i 's/\t$//g' $inputbed

#reduce deletions to consider by intersecting with exons
#note this DOES NOT retain information about the exon being overlapped
bedtools intersect -wa  -a ${inputbed} -b /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${chr}.bed | uniq > chr${chr}.${bin}.del.exonic.bed

#find reciprocal overlaps of more than 80%
# then remove duplicates
# then collapse
bedtools intersect -a chr${chr}.${bin}.del.exonic.bed -b chr${chr}.${bin}.del.exonic.bed -f 0.8 -r -wa -wb | \
	awk -v OFS='\t' 'if($2 != $43 && $3!=$44){print $1,$2,$3,$1":"$2"-"$3,$5,$42,$43,$44,$42":"$43"-"$44,$46}' - | sort -k2,2n -k7,7n | uniq | \
	bedtools groupby -i - -grp 1-4 -c 9 -o collapse | awk -v OFS='\t' '$5 !~ p {print $0};{p=$1":"$2"-"$3}' - > chr${chr}.${bin}.exonic.collapse.key

#write non-overlapping genotypes and collapse overlapping genotypes into the final output file

awk -v OFS='\t' '{print $1":"$2"-"$3,$0}' chr${chr}.${bin}.del.exonic.bed > chr${chr}.${bin}.del.exonic.wIDs.bed

while read -r line ; do id=$(echo $line | cut -f 1 -d " " ) ; if ! grep -q $id chr${chr}.${bin}.exonic.collapse.key  ; \
	then echo $line | cut -f 2- -d " " >> chr${chr}.${bin}.del.exonic.collapsed.bed ; \
	else if ! grep -q $id  chr${chr}.${bin}.del.exonic.collapsed.bed ; \
	then grep $id chr${chr}.${bin}.exonic.collapse.key | cut -f 4-5 | tr "," "\t" | tr "\t" "\n" > ${chr}.${bin}.grep.ids ; \
	alt=$(grep $id chr${chr}.${bin}.exonic.collapse.key | cut -f 4-5 | tr "\t" ","); \
	grep -w -f ${chr}.${bin}.grep.ids chr${chr}.${bin}.del.exonic.wIDs.bed | \
	awk -v OFS='\t' -v alt=$(echo $alt) -v chr=$(echo $chr) '$3 < min {min=$3}; min == 0 {min=$3}; $4 > max {max=$4}; max == 0 {max=$4}; $8 != "." {g8 += $8}; $9 != "." {g9 += $9} ; $10 != "." {g10 += $10} ; $11 != "." {g11 += $11} ; $12 != "." {g12 += $12} ; $13 != "." {g13 += $13} ; $14 != "." {g14 += $14} ; $15 != "." {g15 += $15} ; $16 != "." {g16 += $16} ; $17 != "." {g17 += $17} ; $18 != "." {g18 += $18} ; $19 != "." {g19 += $19} ; $20 != "." {g20 += $20} ; $21 != "." {g21 += $21} ; $22 != "." {g22 += $22} ; $23 != "." {g23 += $23} ; $24 != "." {g24 += $24} ; $25 != "." {g25 += $25} ; $26 != "." {g26 += $26} ; $27 != "." {g27 += $27} ; $28 != "." {g28 += $28} ; $29 != "." {g29 += $29} ; $30 != "." {g30 += $30} ; $31 != "." {g31 += $31} ; $32 != "." {g32 += $32} ; $33 != "." {g33 += $33} ; $14 != "." {g34 += $34} ; $35 != "." {g35 += $35} ; $36 != "." {g36 += $36} ; $37 != "." {g37 += $37} ; $38 != "." {g38 += $38} ; $39 != "." {g39 += $39} ; $40 != "." {g40 += $40} ; $41 != "." {g41 += $41}; $42 != "." {g42 += $42} END {print chr, min,max,"collapsedDEL",alt,"QUAL",g8,g9,g10,g11,g12,g13,g14,g15,g16,g17,g18,g19,g20,g21,g22,g23,g24,g25,g26,g27,g28,g29,g30,g31,g32,g33,g34,g35,g36,g37,g38,g39,g40,g41,g42}' - |\
	sed -E -e :1 -e 's/(^|\t)(\t|$)/\1\.\2/;t1' >> chr${chr}.${bin}.del.exonic.collapsed.bed ; fi ; fi ; done < chr${chr}.${bin}.del.exonic.wIDs.bed

sed -i 's/ /\t/g' chr${chr}.${bin}.del.exonic.collapsed.bed 
```

### bin deletions by length to get a summary of number of exonic deletions of each length

```
awk -v OFS='\t' '$4 != "collapsedDEL" && $3-$2 <= 10 {print "less than 10"};$4 != "collapsedDEL" && $3-$2 > 10 && $3-$2 <= 100 {print "10 to 100 bp"}; $4 != "collapsedDEL" && $3-$2 > 100 && $3-$2 <= 1000 {print "101 to 1000 bp"};$4 != "collapsedDEL" && $3-$2 > 1000 && $3-$2 <= 10000 {print "1000 to 10KB"};$4 != "collapsedDEL" && $3-$2 > 10000 && $3-$2 <= 100000 {print "10KB to 100KB"};$4 != "collapsedDEL" && $3-$2 > 100000 {print "100+KB"}' chr10.bin1.del.exonic.collapsed.bed | sort -k1,1n | uniq -c
awk -v OFS='\t' '$3-$2 <= 10 {print "less than 10"}; $3-$2 > 10 && $3-$2 <= 100 {print "10 to 100 bp"}; $3-$2 > 100 && $3-$2 <= 1000 {print "101 to 1000 bp"}; $3-$2 > 1000 && $3-$2 <= 10000 {print "1000 to 10KB"}; $3-$2 > 10000 && $3-$2 <= 100000 {print "10KB to 100KB"}; $3-$2 > 100000 {print "100+KB"}' chr10.bin1.del.exonic.collapsed.bed | sort -k1,1n | uniq -c
```

### Filter down the collapsed deletions to only those that are EXACTLY contained within exon boundaries

Note that the below commands will include the exon information (not just the deletion information like before)
```
ml bedtools2
for c in 01 02 03 04 05 06 07 08 09 10; do for b in bin1 bin2; do 
qbedtools intersect -wa -wb -F 1 -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${c}.bed -b chr${c}.${b}.del.exonic.collapsed.bed > chr${c}.${b}.del.exact.exonic.collapsed.bed
done; done

#number of unique deletions
wc -l *exact.exonic.collapsed.bed

#number of deletions in each size bin
for c in 01 02 03 04 05 06 07 08 09 10; do for b in bin1 bin2; do echo chr ${c} ${b} ; 
awk -v OFS='\t' '$9-$8 <= 10 {print "less than 10"}; $9-$8 > 10 && $9-$8 <= 100 {print "10 to 100 bp"}; $9-$8 > 100 && $9-$8 <= 1000 {print "101 to 1000 bp"}; $9-$8 > 1000 && $9-$8 <= 10000 {print "1000 to 10KB"}; $9-$8 > 10000 && $9-$8 <= 100000 {print "10KB to 100KB"}; $9-$8 > 100000 {print "100+KB"}' chr${c}.${b}.del.exact.exonic.collapsed.bed | sort -k1,1n | uniq -c
done ; done

#number of unique exons
for c in 01 02 03 04 05 06 07 08 09 10; do for b in bin1 bin2; do
echo In chr ${c} ${b} total number of exons with a deletion:
cut -f 4 chr${c}.${b}.del.exact.exonic.collapsed.bed | sort | uniq | wc -l #number of total exons
echo In chr ${c} ${b} number of exons with a single deletion:
cut -f 4 chr${c}.${b}.del.exact.exonic.collapsed.bed | sort | uniq -u | wc -l #number of exons with 1 deletion
done ; done

#creating summary files that has the counts of dels per exon in the different size categories
for c in 01 02 03 04 05 06 07 08 09 10; do for b in bin1 bin2; do echo chr ${c} ${b} ; 
awk -v OFS='\t' '$9-$8 <= 10 {print $4, "less_than_10"}; $9-$8 > 10 && $9-$8 <= 100 {print $4,"10_to_100bp"}; $9-$8 > 100 && $9-$8 <= 1000 {print $4,"101_to_1000bp"}; $9-$8 > 1000 && $9-$8 <= 10000 {print $4,"1000_to_10KB"}' chr${c}.${b}.del.exact.exonic.collapsed.bed | sort | uniq -c | sed 's/^[_[:space:]]*//' | tr " " "\t" > chr${c}.${b}.del.exact.exonic.collapsed.summary.tsv
done ; done

#printing the sizes of deletions where there's one del per exon
for c in 01 02 03 04 05 06 07 08 09 10; do for b in bin1 bin2; do echo chr ${c} ${b} ; 
awk -v OFS='\t' '$1 == 1 {print $3}' chr${c}.${b}.del.exact.exonic.collapsed.summary.tsv  | sort | uniq -c
done; done

#printing the sizes of deletions where there's multiple dels per exon
for c in 01 02 03 04 05 06 07 08 09 10; do for b in bin1 bin2; do echo chr ${c} ${b} ; 
awk -v OFS='\t' '$1 > 1 {print $3}' chr${c}.${b}.del.exact.exonic.collapsed.summary.tsv  | sort | uniq -c
done; done
```

## Estimate dN/dS
### Pull out sequences for each genome
Will need to get the coordinates for each genome? 

GVCF: info field: ASM_Chr= ASM_End, ASM_Start, ASM_Strand
MAF: s Chr Start Stop Strand Length Sequence
.anchors: has gene names for Sb (but only those genes used as anchors)
```
#make gene IDs list
cut -f 4 ref_Sb313.cds.bed | cut -f 2 -d ";" | sed 's/Parent=//g' | uniq > ref_Sb313.cds.genelist.txt

for i in *.anchors ; do grep -f ref_Sb313.cds.genelist.txt $i > ${i%_anchorwave.anchors}.refGenesAsAnchors.txt ; done

#to figure out how many ref genes show up at least 1x as an anchor:
cut -f 8 *.refGenesAsAnchors.txt | sort | uniq -c | wc -l

#to figure out how many genes show up 1x vs 2x as an anchor:
cut -f 8 *.refGenesAsAnchors.txt | sort | uniq -c | awk '{print $1}' | sort | uniq -c

for i in *refGenesAsAnchors.txt ; do
	echo ${i%.refGenesAsAnchors.txt} had this many ref genes show up at least 1x as an anchor >> refGenesAsAnchors.summary.txt
	cut -f 8 ${i} |  sort | uniq -c | wc -l >> refGenesAsAnchors.summary.txt
	echo ${i%.refGenesAsAnchors.txt} had this many genes show up 1x vs 2x as an anchor >> refGenesAsAnchors.summary.txt 
	cut -f 8 ${i} |  sort | uniq -c | awk '{print $1}' | sort | uniq -c >> refGenesAsAnchors.summary.txt 
done
```
What this tells us is that ~8200 genes out of the 12K are used as anchors for any given genome
Most of those that are used as anchors are only used once

It may help to have a list of genes that are COMPLETELY missing in each genome

```
ml bedtools2
for c in 01 02 03 04 05 06 07 08 09 10 ; do 
	for b in bin1 bin2 ; do
		bedtools intersect -wa -wb -a ../../../ref_Sb313.cds.${c}.bed -b chr${c}.${b}.del.exonic.collapsed.bed -f 1.0 > chr${c}.${b}.missing.exons.collapsed.bed
	done
done	

#colnames: SbCHROM SbStart SbEnd ID QUAL Strand delChrom delStart delStop REF ALT QUAL TdFL TdKS ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25

#1. split alt allele
#2. make sure allele is a deletion
#3. create a list for each genome and add to it when that genome has the deletion allele (denote subgenome)
#Goal: Subgenome, CDS_ID

#13 is TdFL
#so long as the TdFL is not REF or NA, 
#split the alt allele 
#if the genotype is alt 1 & the length of the deletion (length of Ref - length ALT) is greater than/equal to the length of the CDS (End-Start)

awk -v OFS='\t' '$13 != 0 && $13 != "." {split($11,a,/,/); if($13 == 1 && length($10)-length(a[1]) >= $3-$2 ) print $1,$2,$3,$4,$6,"bin1", $7,$8,$9 ; if($13 == 2 && length($10)-length(a[2]) >=$3-$2) print $1,$2,$3,$4,$6,"bin1",$7,$8,$9}'

#creating slurm scripts like: 
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --time=2:00:00
#SBATCH --job-name=TdFL_missing
#SBATCH --output=nova-%x.%j.out
#SBATCH --error=nova-%x.%j.err

for i in *bin1.missing.exons.collapsed.bed ; do awk -v OFS='\t' '$13 != 0 && $13 != "." {split($11,a,/,/); if($13 == 1 && length($10)-length(a[1]) >= $3-$2 ) print $1,$2,$3,$4,$6,"bin1", $7,$8,$9 ; if($13 == 2 && length($10)-length(a[2]) >=$3-$2) print $1,$2,$3,$4,$6,"bin1",$7,$8,$9}' ${i} >> TdFL.missing.exons.collapsed.bed
done
for i in *bin2.missing.exons.collapsed.bed ; do awk -v OFS='\t' '$13 != 0 && $13 != "." {split($11,a,/,/); if($13 == 1 && length($10)-length(a[1]) >= $3-$2 ) print $1,$2,$3,$4,$6,"bin2", $7,$8,$9 ; if($13 == 2 && length($10)-length(a[2]) >=$3-$2) print $1,$2,$3,$4,$6,"bin2",$7,$8,$9}' ${i} >> TdFL.missing.exons.collapsed.bed
done

###
#To make copies of the script for other genomes:
for g in TdKS ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
	cp slurm_createMissingExonByGenome.TdFL.sh slurm_createMissingExonByGenome.${g}.sh
	sed -i "s/TdFL/${g}/g" slurm_createMissingExonByGenome.${g}.sh
done
sed -i 's/13/14/g' slurm_createMissingExonByGenome.TdKS.sh
sed -i 's/13/15/g' slurm_createMissingExonByGenome.ZdGigi_4to1.sh 
sed -i 's/13/16/g' slurm_createMissingExonByGenome.ZdMomo_4to1.sh 
sed -i 's/13/17/g' slurm_createMissingExonByGenome.ZhRIMHU001.sh 
sed -i 's/13/18/g' slurm_createMissingExonByGenome.ZmB73.sh 
sed -i 's/13/19/g' slurm_createMissingExonByGenome.ZmB97.sh 
sed -i 's/13/20/g' slurm_createMissingExonByGenome.ZmCML103.sh 
sed -i 's/13/21/g' slurm_createMissingExonByGenome.ZmCML228.sh
sed -i 's/13/22/g' slurm_createMissingExonByGenome.ZmCML247.sh 
sed -i 's/13/23/g' slurm_createMissingExonByGenome.ZmCML277.sh 
sed -i 's/13/24/g' slurm_createMissingExonByGenome.ZmCML322.sh 
sed -i 's/13/25/g' slurm_createMissingExonByGenome.ZmCML333.sh 
sed -i 's/13/26/g' slurm_createMissingExonByGenome.ZmCML52.sh 
sed -i 's/13/27/g' slurm_createMissingExonByGenome.ZmCML69.sh 
sed -i 's/13/28/g' slurm_createMissingExonByGenome.ZmHP301.sh 
sed -i 's/13/29/g' slurm_createMissingExonByGenome.ZmIL14H.sh 
sed -i 's/13/30/g' slurm_createMissingExonByGenome.ZmKi11.sh 
sed -i 's/13/31/g' slurm_createMissingExonByGenome.ZmKi3.sh 
sed -i 's/13/32/g' slurm_createMissingExonByGenome.ZmKy21.sh 
sed -i 's/13/33/g' slurm_createMissingExonByGenome.ZmM162W.sh 
sed -i 's/13/34/g' slurm_createMissingExonByGenome.ZmM37W.sh 
sed -i 's/13/35/g' slurm_createMissingExonByGenome.ZmMS71.sh 
sed -i 's/13/36/g' slurm_createMissingExonByGenome.ZmMo18W.sh 
sed -i 's/13/37/g' slurm_createMissingExonByGenome.ZmNC350.sh 
sed -i 's/13/38/g' slurm_createMissingExonByGenome.ZmNC358.sh 
sed -i 's/13/39/g' slurm_createMissingExonByGenome.ZmOh43.sh 
sed -i 's/13/40/g' slurm_createMissingExonByGenome.ZmOh7b.sh 
sed -i 's/13/41/g' slurm_createMissingExonByGenome.ZmP39.sh 
sed -i 's/13/42/g' slurm_createMissingExonByGenome.ZmTx303.sh 
sed -i 's/13/43/g' slurm_createMissingExonByGenome.ZmTzi8.sh 
sed -i 's/13/44/g' slurm_createMissingExonByGenome.ZvTIL01.sh 
sed -i 's/13/45/g' slurm_createMissingExonByGenome.ZvTIL11.sh 
sed -i 's/13/46/g' slurm_createMissingExonByGenome.ZxTIL18.sh 
sed -i 's/13/47/g' slurm_createMissingExonByGenome.ZxTIL25.sh

#So to get the CDS_ID list, you'd just need to do:
cut -f 4 TdFL.missing.exons.collapsed.bed | tr ";" "\t" | cut -f 1 | sed 's/ID=//g'  
#There could be multiple deletion alleles that cover the exon, so may not just CDS ID may not just show up 1 or 2

cut -f 4,6 TdFL.missing.exons.collapsed.bed | tr ";" "\t" | cut -f 1,4 | sed 's/ID=//g'   | sort | uniq -c | awk '{print $1}' | sort | uniq -c 
```
This is great, but it would be helpful to have a file that had the number of exons in each gene model
Maybe we to get at extent of fractionation, we say something about the % of basepairs that intersect with a deletion /gene model?



Let's try to subset out the genes of interest using the `maf_split` tools from `phast`
And let's do it one gene at a time

Step 0: create a file with the unique Gene ID names
The first field is the number of exons in that gene model
The second field is the name of the gene model
```
cut -f 4 /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed | cut -f 2 -d ";" | sed "s/Parent=//g" | sort | uniq -c |sed 's/^[_[:space:]]*//; s/[_[:space:]]*$//' | tr " " "\t" > ref_Sb313.geneids.txt
```

Create a directory for each gene
```
mkdir gene_alignments
cd gene_alignments

mkdir Sobic.001G000200.1.v3.1
cd Sobic.001G000200.1.v3.1
grep -w Sobic.001G000200.1.v3.1 /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed | awk -v OFS='\t' '{print "Chr"$1,$2,$3,$4,1,$6}' - > Sobic.001G000200.1.v3.1.bed

awk -v OFS='\t' '!i++ {min = $2; max = $3} { min = (min < $2) ? min : $2 ; max = (max > $3) ? max : $3} END {print $1,min,max,$4,$5,$6}' Sobic.001G000200.1.v3.1.bed > Sobic.001G000200.1.v3.1.gene.bed 

conda activate phast
while read line ; do 
	echo ${line} | tr " " "\t" > temp.bed
	for g in TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
		maf_parse --features temp.bed ../../${g}_Chr01.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Sobic.001G000200.1.v3.1.bin1.${g}.maf
		maf_parse --features temp.bed ../../${g}_Chr01.bin2.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Sobic.001G000200.1.v3.1.bin2.${g}.maf
		if [ ${g} == TdFL ] ; then 
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 ~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin1.${g}.maf |tr "\t" "\n" >> Sb313.Sobic.001G000200.1.v3.1.bin1.cds.fasta
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s"  && $2 ~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin2.${g}.maf |tr "\t" "\n" >> Sb313.Sobic.001G000200.1.v3.1.bin2.cds.fasta
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin1.${g}.maf |tr "\t" "\n" >> ${g}.Sobic.001G000200.1.v3.1.bin1.cds.fasta
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s"  && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin2.${g}.maf |tr "\t" "\n" >> ${g}.Sobic.001G000200.1.v3.1.bin2.cds.fasta
			else
				awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin1.${g}.maf |tr "\t" "\n" >> ${g}.Sobic.001G000200.1.v3.1.bin1.cds.fasta
				awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin2.${g}.maf |tr "\t" "\n" >> ${g}.Sobic.001G000200.1.v3.1.bin2.cds.fasta
		fi
	done
done < Sobic.001G000200.1.v3.1.bed

for i in {1..11} ; do grep -A 1 "CDS.${i};" TdFL.Sobic.001G000200.1.v3.1.bin1.cds.fasta | awk -v OFS='\t' '$1 !~ "ID=" {print $0}' >> temp.TdFL.Sobic.001G000200.1.v3.1.bin1.cds.fasta ; done
echo \>Sobic.001G000200.1.v3.1\;TdFL\;+ >> Sobic.001G000200.1.v3.1.bin1.cds.fasta
tr -d "\n\r" < temp.TdFL.Sobic.001G000200.1.v3.1.bin1.cds.fasta >> Sobic.001G000200.1.v3.1.bin1.cds.fasta

for g in Sb313 TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25; do 
	for i in {1..11}; do 
		grep -A 1 "CDS.${i};" ${g}.Sobic.001G000200.1.v3.1.bin1.cds.fasta | awk -v OFS='\t' '$1 !~ "ID=" {print $0}' >> temp.${g}.Sobic.001G000200.1.v3.1.bin1.cds.fasta ; 
	done
	printf '\n%s\n' \>Sobic.001G000200.1.v3.1\;${g}\;+ >> Sobic.001G000200.1.v3.1.bin1.cds.fasta
	tr -d "\n\r" < temp.${g}.Sobic.001G000200.1.v3.1.bin1.cds.fasta >> Sobic.001G000200.1.v3.1.bin1.cds.fasta
done

muscle -in Sobic.001G000200.1.v3.1.bin1.cds.fasta  -out Sobic.001G000200.1.v3.1.bin1.cds.aln.fasta  -maxiters 2
muscle -in Sobic.001G000200.1.v3.1.bin2.cds.fasta  -out Sobic.001G000200.1.v3.1.bin2.cds.aln.fasta  -maxiters 2 #throws an error because there's no alignment for bin2 for this gene


for g in TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
		maf_parse --features Sobic.001G000200.1.v3.1.gene.bed  ../../${g}_Chr01.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Sobic.001G000200.1.v3.1.bin1.${g}.gene.maf
		maf_parse --features Sobic.001G000200.1.v3.1.gene.bed  ../../${g}_Chr01.bin2.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Sobic.001G000200.1.v3.1.bin2.${g}.gene.maf
		if [ ${g} == TdFL ] ; then 
			awk -v id=$(echo Sobic.001G000200.1.v3.1 ) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin1.${g}.gene.maf |tr "\t" "\n" >> Sobic.001G000200.1.v3.1.bin1.gene.fasta
			awk -v id=$(echo Sobic.001G000200.1.v3.1 ) -v OFS='\t' '$1=="s"  {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin2.${g}.gene.maf |tr "\t" "\n" >> Sobic.001G000200.1.v3.1.bin2.gene.fasta
			else
				awk -v id=$(echo Sobic.001G000200.1.v3.1 ) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin1.${g}.gene.maf |tr "\t" "\n" >> Sobic.001G000200.1.v3.1.bin1.gene.fasta
				awk -v id=$(echo Sobic.001G000200.1.v3.1 ) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Sobic.001G000200.1.v3.1.bin2.${g}.gene.maf |tr "\t" "\n" >> Sobic.001G000200.1.v3.1.bin2.gene.fasta
		fi
	done
rm temp.*

ml muscle
muscle -in Sobic.001G000200.1.v3.1.bin1.gene.fasta  -out Sobic.001G000200.1.v3.1.bin1.gene.aln.fasta  -maxiters 2
muscle -in Sobic.001G000200.1.v3.1.bin2.gene.fasta  -out Sobic.001G000200.1.v3.1.bin2.gene.aln.fasta  -maxiters 2 #throws an error because there's no alignment for bin2 for this gene

conda deactivate
```
Do I need to split out to concatenate the exons into CDS for each genome... yes
But that doesn't need to happen for the gene maf parse steps

Now to make it a script that can be looped through by 12K gene IDs

```
#!/bin/bash

geneid=$1
exoncnt=$2

#make the directory for the gene id and move there
cd /work/LAS/mhufford-lab/snodgras/Fractionation/AnchorWave_output/tripsacinae-sb_split_mafs/gene_alignments
mkdir ${geneid}
cd ${geneid}

echo Directory made...

#create bed files for CDS and gene
grep -w ${geneid} /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.bed | awk -v OFS='\t' '{print "Chr"$1,$2,$3,$4,1,$6}' - > ${geneid}.bed
awk -v OFS='\t' '!i++ {min = $2; max = $3} { min = (min < $2) ? min : $2 ; max = (max > $3) ? max : $3} END {print $1,min,max,$4,$5,$6}' ${geneid}.bed > ${geneid}.gene.bed 

echo Bed files made...

#pull out sequences from MAFs
## CDS
conda activate phast
while read line ; do 
	echo ${line} | tr " " "\t" > temp.bed
	CHR=$(echo $line | cut -f 1 -d " ")
	for g in TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
		maf_parse --features temp.bed ../../${g}_${CHR}.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.${geneid}.v3.1.bin1.${g}.maf
		maf_parse --features temp.bed ../../${g}_${CHR}.bin2.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.${geneid}.bin2.${g}.maf
		if [ ${g} == TdFL ] ; then 
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 ~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin1.${g}.maf |tr "\t" "\n" >> Sb313.${geneid}.bin1.cds.fasta
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s"  && $2 ~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin2.${g}.maf |tr "\t" "\n" >> Sb313.${geneid}.bin2.cds.fasta
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin1.${g}.maf |tr "\t" "\n" >> ${g}.${geneid}.bin1.cds.fasta
			awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s"  && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin2.${g}.maf |tr "\t" "\n" >> ${g}.${geneid}.bin2.cds.fasta
			else
				awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin1.${g}.maf |tr "\t" "\n" >> ${g}.${geneid}.bin1.cds.fasta
				awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin2.${g}.maf |tr "\t" "\n" >> ${g}.${geneid}.bin2.cds.fasta
		fi
	done
done < ${geneid}.bed

echo CDS sequences pulled out of MAFs...

## gene
for g in TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
		CHR=$(head -n 1 ${geneid}.gene.bed | cut -f 1 -d " ")
		maf_parse --features ${geneid}.gene.bed  ../../${g}_${CHR}.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.${geneid}.v3.1.bin1.${g}.gene.maf
		maf_parse --features ${geneid}.gene.bed  ../../${g}_${CHR}.bin2.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.${geneid}.bin2.${g}.gene.maf
		if [ ${g} == TdFL ] ; then 
			awk -v id=$(echo ${geneid} ) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin1.${g}.gene.maf |tr "\t" "\n" >> ${geneid}.bin1.gene.fasta
			awk -v id=$(echo ${geneid} ) -v OFS='\t' '$1=="s"  {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.${g}.gene.maf |tr "\t" "\n" >> ${geneid}.bin2.gene.fasta
			else
				awk -v id=$(echo ${geneid} ) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin1.${g}.gene.maf |tr "\t" "\n" >> ${geneid}.bin1.gene.fasta
				awk -v id=$(echo ${geneid} ) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.${geneid}.bin2.${g}.gene.maf |tr "\t" "\n" >> ${geneid}.bin2.gene.fasta
		fi
done

echo Gene sequences pulled out of MAFs...

#create combined CDS fasta
for g in Sb313 TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25; do 
	for (( i=1; i<=${exonct}; c++)); do 
		grep -A 1 "CDS.${i};" ${g}.${geneid}.bin1.cds.fasta | awk -v OFS='\t' '$1 !~ "ID=" {print $0}' >> temp.${g}.${geneid}.bin1.cds.fasta ;
		grep -A 1 "CDS.${i};" ${g}.${geneid}.bin2.cds.fasta | awk -v OFS='\t' '$1 !~ "ID=" {print $0}' >> temp.${g}.${geneid}.bin2.cds.fasta ; 
	done
	printf '\n%s\n' \>${geneid}\;${g}\;+ >> ${geneid}.bin1.cds.fasta
	tr -d "\n\r" < temp.${g}.${geneid}.bin1.cds.fasta >> ${geneid}.bin1.cds.fasta
	printf '\n%s\n' \>${geneid}\;${g}\;+ >> ${geneid}.bin2.cds.fasta
	tr -d "\n\r" < temp.${g}.${geneid}.bin1.cds.fasta >> ${geneid}.bin2.cds.fasta
done

echo CDS sequences combined into fasta...

#remove temporary files
rm temp*

conda deactivate

#create muscle alignments
#throws an error if there's no alignment
ml muscle
muscle -in ${geneid}.bin1.gene.fasta  -out ${geneid}.bin1.gene.aln.fasta  -maxiters 2
muscle -in ${geneid}.bin2.gene.fasta  -out ${geneid}.bin2.gene.aln.fasta  -maxiters 2 
muscle -in ${geneid}.bin1.cds.fasta  -out ${geneid}.bin1.cds.aln.fasta  -maxiters 2
muscle -in ${geneid}.bin2.cds.fasta  -out ${geneid}.bin2.cds.aln.fasta  -maxiters 2 

```


Need to do a check to see if I should push it through to the dN/dS step
See if there's a bunch of missing introns

```
ml bioawk
bioawk -c fastx '{ print $name, length($seq) }' < sequences.fa

```

Could skip dN/dS calculations and use what's already been published by Yin et al 2022 MBE
```
#get the supplemental table 3
#copy the sorghum syntelog ids into a text file called "Yin2022MBE.SbGeneIDs.txt"

grep -f Yin2022MBE.SbGeneIDs.txt ref_Sb313.geneids.txt > shared_geneIDs.txt

```
Of out 12,169 ref gene models and the 4,578 genes in Yin et al 2022, 3,195 genes are shared
 





This is how you do it, but with disregard to how many exons you're pulling out (aka pulling all exons out and putting them all in the same fasta)
```
ml singularity
singularity run https://depot.galaxyproject.org/singularity/phast:1.5--hec16e2b_5

head -n 1 ../ref_Sb313.cds.bed | awk -v OFS='\t' '{print "Chr"$1,$2,$3,$4,1,$6}' - > test.bed

maf_parse --features test.bed TdFL_Chr01.bin1.maf #has to be one chromosome at a time
sed -i 's/Chr01/Chr01.Sb313/g' test.TdFL_Chr01.bin1.exon1.maf 
sed -i 's/chr1/chr1.TdFL/g' test.TdFL_Chr01.bin1.exon1.maf 

maf_parse --features ../test.bed TdKS_Chr01.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed 's/chr/TdKS.chr/g' > test.TdKS_Chr01.bin1.exon1.maf 
awk -v id=$(cut -f 4 ../test.bed) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' test.TdKS_Chr01.bin1.exon1.maf | tr "\t" "\n" > test.TdKS_Chr01.bin1.exon1.fa

head -n 3 ../ref_Sb313.cds.bed | awk -v OFS='\t' '{print "Chr"$1,$2,$3,$4,1,$6}' - > test.bed

maf_parse --features ../test.bed TdFL_Chr01.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed 's/chr/TdFL.chr/g' > test.TdFL_Chr01.bin1.exon1.maf 

while read line ; do 
	echo $line | tr " " "\t" > temp.bed
	maf_parse --features temp.bed TdKS_Chr01.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed 's/chr/TdKS.chr/g' > temp.maf
	awk -v id=$(cut -f 4 temp.bed) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.maf | tr "\t" "\n" >> test.TdKS_Chr01.bin1.fa
	rm temp*
done < ../test.bed

ml muscle
muscle -in Chr01.bin1.fasta -out Chr01.bin1.aln.fasta -maxiters 2

```
Let's turn it into a script that can run on each chromosome/bin

```
conda activate phast

while read line ; do 
	echo $line | tr " " "\t" |awk -v OFS='\t' '{print "Chr"$1,$2,$3,$4,1,$6}' - > temp.Chr01.bed
	for g in TdFL TdKS ZnPI615697_4to1 ZdGigi_4to1 ZdMomo_4to1 ZhRIMHU001 ZmB73 ZmB97 ZmCML103 ZmCML228 ZmCML247 ZmCML277 ZmCML322 ZmCML333 ZmCML52 ZmCML69 ZmHP301 ZmIL14H ZmKi11 ZmKi3 ZmKy21 ZmM162W ZmM37W ZmMS71 ZmMo18W ZmNC350 ZmNC358 ZmOh43 ZmOh7b ZmP39 ZmTx303 ZmTzi8 ZvTIL01 ZvTIL11 ZxTIL18 ZxTIL25 ; do
		maf_parse --features temp.Chr01.bed ${g}_Chr01.bin1.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Chr01.bin1.${g}.maf
		maf_parse --features temp.Chr01.bed ${g}_Chr01.bin2.maf | sed 's/Chr/Sb313.Chr/g' | sed "s/chr/${g}.chr/g" > temp.Chr01.bin2.${g}.maf
		if [ ${g} == TdFL ] ; then 
			awk -v id=$(cut -f 4 temp.Chr01.bed) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.Chr01.bin1.${g}.maf |tr "\t" "\n" >> Chr01.bin1.fasta
			awk -v id=$(cut -f 4 temp.Chr01.bed) -v OFS='\t' '$1=="s" {print ">"id";"$2";"$3";"$5, $7}' temp.Chr01.bin2.${g}.maf |tr "\t" "\n" >> Chr01.bin2.fasta
			else
				awk -v id=$(cut -f 4 temp.Chr01.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Chr01.bin1.${g}.maf |tr "\t" "\n" >> Chr01.bin1.fasta
				awk -v id=$(cut -f 4 temp.Chr01.bed) -v OFS='\t' '$1=="s" && $2 !~ /Sb313/ {print ">"id";"$2";"$3";"$5, $7}' temp.Chr01.bin2.${g}.maf |tr "\t" "\n" >> Chr01.bin2.fasta
		fi
	done
	rm temp.Chr01*
done < /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.01.bed

conda deactivate
           
        
for i in {2..10} ; do 
cp slurm_11.1.maf2fasta.sh slurm_11.${i}.maf2fasta.sh 
if [ i != 10 ] ; then sed -i "s/Chr01/Chr0${i}/g" slurm_11.${i}.maf2fasta.sh ; else sed -i "s/Chr01/Chr${i}/g" slurm_11.${i}.maf2fasta.sh ; fi
if [ i != 10 ] ; then sed -i "s/cds\.01\.bed/cds\.0${i}\.bed" slurm_11.${i}.maf2fasta.sh ; else sed -i "s/cds\.01\.bed/cds\.${i}\.bed" slurm_11.${i}.maf2fasta.sh ; fi
done
 ```
            
 










_OLD METHOD_
### test
```
mkdir test_vcf2bed
grep "^##" -v chr10.bin2.indelonly_reformatted.vcf > test_vcf2bed/chr10.bin2.headerless.indelonly.vcf

cd test_vcf2bed

#this takes the headerless vcf and checks the alt allele; at least 1 ALT has to be smaller in length than the REF allele
head -n 1 chr10.bin2.headerless.indelonly.vcf > chr10.bin2.delonly.vcf
cat chr10.bin2.headerless.indelonly.vcf | \
awk -v OFS='\t' '{split($5,a,/,/); if(length(a[1]) < length($4) || (length(a[2]) < length($4) && a[2] != "")) print $0}' - >> chr10.bin2.delonly.vcf

#but sometimes this can be insertions+del or deletions of different lengths
#so next need each ALT on its own line and change the genotypes accordingly? 

#makes it a bed with the stop being the largest possible of the two alt alleles
awk -v OFS='\t' '{split($5,a,/,/); 
	if(length($4) - length(a[1]) >= length($4) - length(a[2])) 
		s=(length($4) - length(a[1])); else s=length($4) - length(a[2]);
		print $1,$2,$2+s,$4,$5,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44}' chr10.bin2.delonly.vcf > chr10.bin2.delonly.bed

#how may instances of insertions and deletions as alt alleles for the same ref pos
awk -v OFS='\t' '{split($5,a,/,/); if(length(a[1]) < length($4) && length(a[2]) < length($4)) print $0}' chr10.bin2.delonly.bed | wc -l 
#number that are only deletions for both ALT alleles or only have 1 deletion ALT allele
258342
wc -l chr10.bin2.delonly.bed #total number of ref positions
277868 #includes header line

19525 instances where there's an insertion and deletion alt allele
```


## Compare numbers of deletions that are
*Test this on a test VCF first*
### exonic vs. non-exonic
#### test
```
ml bedtools2

bedtools intersect -wa -wb -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.10.bed -b chr10.bin2.delonly.bed > chr10.bin2.del.exonic.bed
bedtools intersect -v -wa -b /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.10.bed -a chr10.bin2.delonly.bed > chr10.bin2.del.nonexonic.bed

#how many non-exonic deletions
wc -l chr10.bin2.del.nonexonic.bed #270620

#how many exonic deletions

cut -f 7-47 chr10.bin2.del.exonic.bed | sort -k2,3 | uniq | wc -l 
#7247

```

### single deletions (biallelic) vs. nested (multi-allelic)
#### test
```
#how many deletions per exon?
cut -f 1-6 chr10.bin2.del.exonic.bed | sort | uniq | wc -l  
#3825 unique exons have a deletion
wc -l chr10.bin2.del.exonic.bed 
# 26893 total number of exons x deletions intersecting them

cut -f 4 chr10.bin2.del.exonic.bed | sort | uniq -u | wc -l #should be the number of exons that only intersect with 1 deletion (327)
#get the number of times X number of deletions overlap within a single exon
echo Number_of_Exons Number_of_OverLapping_Deletions > test.exonicDelOverlap.txt
cut -f 4 chr10.bin2.del.exonic.bed | sort | uniq -c | sed 's/ //g' | sed 's/ID/ ID/g' | cut -f 1 -d " " | sort | uniq -c  >>test.exonicDelOverlap.txt
#So 3498 exons have multiple deletion overlaps

### NOTE: -c option in bedtools intersect will count the number of B entries that overlap with A entries
```

### within exon boundaries vs. extending past
```
ml bedtools2

bedtools intersect -wa -wb -F 1 -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.10.bed -b chr10.bin2.delonly.bed > chr10.bin2.del.exact.exonic.bed
wc -l chr10.bin2.del.exact.exonic.bed 
#4505 chr10.bin2.del.exact.exonic.bed #total number of exons x deletions intersecting within boundaries
cut -f 4 chr10.bin2.del.exact.exonic.bed | sort | uniq | wc -l #655 exons that have deletions within boundaries of start and stop
cut -f 4 chr10.bin2.del.exact.exonic.bed | sort | uniq -u | wc -l #161 exons are within boundaries of an exon and are the only deletion within that exon
cut -f 4 chr10.bin2.del.exact.exonic.bed  | sort | uniq -c | sed 's/ //g' | sed 's/ID/ ID/g' | cut -f 1 -d " " | sort | uniq -c

```

### in-frame deletions vs. frame shift
```
#Need to figure out how to get the length of the deletion here... might be easier to do $9-$8
#awk -v OFS='\t' '{if($10 % 3 == 0)print $0}' chr10.bin2.del.exonic.bed > chr10.bin2.del.exonic.divisibleBy3.bed

#Does the transition from vcf to bed cause an issue for 0 or 1 indexing? 
#How will that affect the length calculation?
awk -v OFS='\t' '{l = $9-$8 ; if(l % 3 == 0)print $0}' chr10.bin2.del.exact.exonic.bed > chr10.bin2.del.exact.exonic.divisibleBy3.bed

#to see how many unique exons have multiple in frame deletions
cut -f 4 chr10.bin2.del.exact.exonic.divisibleBy3.bed | sort | uniq -c #ranges from 1-36, 456 unique exons

#to find the frame shifts
awk -v OFS='\t' '{l = $9-$8 ; if(l % 3 != 0)print $0}' chr10.bin2.del.exact.exonic.bed > chr10.bin2.del.exact.exonic.frameshift.bed
cut -f 4 chr10.bin2.del.exact.exonic.frameshift.bed | sort | uniq -c #1-50something, 525 unique exons

#How many exons are in both lists (have both in frame and frame shift deletions)
cut -f 4 chr10.bin2.del.exact.exonic.divisibleBy3.bed | grep -f - chr10.bin2.del.exact.exonic.frameshift.bed | cut -f 4 | sort | uniq | wc -l
#326 (between 60 and 70% of the exons on each list)


```

### by exon order (more deletions towards the end of the gene model)

## Summarize the above code in unix so smaller files can be uploaded to R to avoid slowing

```10.vcfFilteringAndReformmating.sh
#!/bin/bash
vcf=$1 #like: chr10.bin2.headerless.indelonly.vcf 
output=${vcf%.headerless.indelonly.vcf}
chr=$2

#filter out insertions
head -n 1 ${vcf} > ${output}.delonly.vcf
cat ${vcf}| \
awk -v OFS='\t' '{split($5,a,/,/); if(length(a[1]) < length($4) || (length(a[2]) < length($4) && a[2] != "")) print $0}' - >> ${output}.delonly.vcf

#makes it a bed with the stop being the end of the REF allele
awk -v OFS='\t' '{split($5,a,/,/); 
		print $1,$2,$2+length($4),$4,$5,$6,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45}' ${output}.delonly.vcf > ${output}.delonly.bed

ml bedtools2

#create exonic 
bedtools intersect -wa -wb -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${chr}.bed -b ${output}.delonly.bed > ${output}.del.exonic.bed

#create exact exonic
bedtools intersect -wa -wb -F 1 -a /work/LAS/mhufford-lab/snodgras/Fractionation/ref_Sb313.cds.${chr}.bed -b ${output}.delonly.bed > ${output}.del.exact.exonic.bed
```

To get the values:
```
#number exonic dels
for i in *del.exonic.bed ; do echo $i ; cut -f 7-48 $i | sort -k2,3 | uniq | wc -l ; done
#number of exons with deletion
for i in *del.exonic.bed ; do echo $i ; cut -f 1-6 $i | sort |uniq | wc -l ; done
#number of single deletions
for i in *del.exonic.bed ; do echo $i ; cut -f 4 $i| sort | uniq -u | wc -l ; done
#number of overlapping deletions
for i in *del.exonic.bed ; do cut -f 4 $i | sort | uniq -c | sed 's/ //g' | sed 's/ID/ ID/g' | cut -f 1 -d " " | sort | uniq -c  >> ${i%.bed}.overlap.txt ;done
awk '{print $1}' *del.exonic.overlap.txt #omit first line and sum the rest
#number of deletions within the boundaries of an exon
for i in *del.exact.exonic.bed ; do wc -l $i ; done

#number in frame deletions
for i in *exact.exonic.bed; do awk -v OFS='\t' '{l = $9-$8 ; if(l % 3 == 0)print $0}' $i > ${i%.bed}.divisibleby3.bed ;done
wc -l *divisibleby3.bed 

#number of frameshift deletions
for i in *exact.exonic.bed; do awk -v OFS='\t' '{l = $9-$8 ; if(l % 3 != 0)print $0}' $i > ${i%.bed}.frameshift.bed ;done
wc -l *frameshift.bed 
```

## Add deletion type to the fractionation status matrix











For looking at ref gene models in R

```
15111 Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt
cut -f 1 Av_Sb313_ExactMatchOnly.OverlapWithCuratedSet.txt | cut -f 1-2 -d "." | sort | uniq | wc -l 
   12169 #Just uniq gene models?
   
#How to deal with isoforms:
#keep the one with the most exons
#if a tie, pick one at random
#This will be easiest to do in R
 
```
```
#To get the number of deletions
for j in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10; do echo $j ; for i in ZmB73 ZmB97 ZmHP301 ZmIL14H ZmKy21 ZmM162W ZmM37W ZmMo18W ZmMS71 ZmOh43 ZmOh7b ZmP39 ZmTx303; do wc -l Sb313_${i}_${j}.bin2.dels.bed | cut -f 1 -d " " ; done ;done

#To get the sum of the deletion lengths
for j in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 ; 
do echo $j ; for i in TdFL	ZdGigi	ZdMomo	ZhRIMHU001	ZnPI615697	ZvTIL01	ZvTIL11	ZxTIL18	ZxTIL25	ZmCML103	ZmCML228	ZmCML247	ZmCML277	ZmCML322	ZmCML333	ZmCML52	ZmCML69	ZmKi11	ZmKi3	ZmNC350	ZmNC358	ZmTzi8	ZmB73	ZmB97	ZmHP301	ZmIL14H	ZmKy21	ZmM162W	ZmM37W	ZmMo18W	ZmMS71	ZmOh43	ZmOh7b	ZmP39	ZmTx303 ; 
do awk '{sum += $8 } END {print sum}' Sb313_${i}_${j}.bin1.dels.bed ; done ; done

for j in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 ; 
do echo $j ; for i in TdFL	ZdGigi	ZdMomo	ZhRIMHU001	ZnPI615697	ZvTIL01	ZvTIL11	ZxTIL18	ZxTIL25	ZmCML103	ZmCML228	ZmCML247	ZmCML277	ZmCML322	ZmCML333	ZmCML52	ZmCML69	ZmKi11	ZmKi3	ZmNC350	ZmNC358	ZmTzi8	ZmB73	ZmB97	ZmHP301	ZmIL14H	ZmKy21	ZmM162W	ZmM37W	ZmMo18W	ZmMS71	ZmOh43	ZmOh7b	ZmP39	ZmTx303 ; 
do awk '{sum += $8 } END {print sum}' Sb313_${i}_${j}.bin2.dels.bed ; done ; done

#To get the sum of the bps in the GVCFs
for j in Chr01 ; do for i in TdFL	ZdGigi	ZdMomo	ZhRIMHU001	ZnPI615697	ZvTIL01	ZvTIL11	ZxTIL18	ZxTIL25	ZmCML103	ZmCML228	ZmCML247	ZmCML277	ZmCML322	ZmCML333	ZmCML52	ZmCML69	ZmKi11	ZmKi3	ZmNC350	ZmNC358	ZmTzi8	ZmB73	ZmB97	ZmHP301	ZmIL14H	ZmKy21	ZmM162W	ZmM37W	ZmMo18W	ZmMS71	ZmOh43	ZmOh7b	ZmP39	ZmTx303 ; 
do zcat Sb313_${i}_${j}.bin1.gvcf.gz | grep -v "^#" | awk -v OFS='\t' '{print $2,$2+length($4)}' | awk -v OFS='\t' '{sum += $2-$1} END {print sum}'
done
done
```
