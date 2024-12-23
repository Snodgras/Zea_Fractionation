#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=12:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

for i in split_mafs_4genomeTrial/unzipped_gvcf/*forEitherBin.exonsNotInGVCF.IDs.txt ; do cat $i >> split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.IDs.temp ;done
sort split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.IDs.temp | uniq > split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.IDs.txt
rm split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.IDs.temp

grep -w -f split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.IDs.txt split_mafs_4genomeTrial/chr10.nodel.forEitherBin.bed > split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.bed

for i in split_mafs_4genomeTrial/unzipped_gvcf/*forEitherBin.exonsInGVCF.IDs.txt ; do cat $i >> split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.IDs.temp ;done
sort split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.IDs.temp | uniq > split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.IDs.txt
rm split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.IDs.temp

grep -w -f split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.IDs.txt split_mafs_4genomeTrial/chr10.nodel.forEitherBin.bed > split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.bed

mkdir QC_exon_fasta

ml bedtools2
for i in assemblies_final/Td-FL assemblies_final/Zd-Gigi assemblies_final/Zv-TIL01 NAM-assemblies/Zm-B73 ; do
	bedtools getfasta -fi ${i}*.fasta -bed split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsNotInGVCF.bed -name -fo QC_exon_fasta/exonsNotInGVCF_${i#*/}.fasta
	bedtools getfasta -fi ${i}*.fasta -bed split_mafs_4genomeTrial/unzipped_gvcf/ALL.nodel.forEitherBin.exonsInGVCF.bed -name -fo QC_exon_fasta/exonsInGVCF_${i#*/}.fasta
done
