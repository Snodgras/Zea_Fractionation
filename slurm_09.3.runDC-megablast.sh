#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --ntasks=36 
#SBATCH --mem=350GB 
#SBATCH --time=12:00:00
#SBATCH --mail-user=snodgras@iastate.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

#!/bin/bash

module load blast-plus

mkdir blast_output

# have all your masked blast databases within the same directory

#for database in blastdb_files/*.nhr
for database in blastdb_files/Td-FL*.nhr blastdb_files/Zd-Gigi*.nhr blastdb_files/Zv-TIL01*.nhr blastdb_files/Zm-B73*.nhr
        do

 n=${database%.*}

blastn -task dc-megablast -outfmt 6 -query QC_exon_fasta/exonsNotInGVCF_Sb313.chr10.fasta -db ${database%.*} -num_threads 4 -out "blast_output/${n#blastdb_files/}_exonsNotInGVCF.txt"
blastn -task dc-megablast -outfmt 6 -query QC_exon_fasta/exonsInGVCF_Sb313.chr10.fasta -db ${database%.*} -num_threads 4 -out "blast_output/${n#blastdb_files/}_exonsInGVCF.txt"

        done

