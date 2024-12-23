#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=12:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

module load blast-plus

mkdir blastdb_files

#looping blastdb:

for sample in assemblies_final/*.fasta
        do
                echo $sample
                describer=$(echo ${sample#assemblies_final/} | cut -f 1-2 -d "-")
                echo $describer

makeblastdb -in ${sample} -dbtype nucl -out blastdb_files/${describer}

done

for sample in NAM-assemblies/*.fasta
        do
                echo $sample
		describer=$(echo ${sample#NAM-assemblies/} | cut -f 1 -d ".")
                echo $describer

makeblastdb -in ${sample} -dbtype nucl -out blastdb_files/${describer}

done
