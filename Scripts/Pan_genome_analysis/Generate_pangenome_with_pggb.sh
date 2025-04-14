#!/bin/bash 

# load pggb in your environment 
source /opt/share/software/scs/appStore/modules/init/profile.sh
module load pggb/v0.5.3


# define variables: 
input_file="/path/to/sloci/.fasta"
flanking_nucleotides="number_of_flanking_bp"
genomes="/path/to/pggb/s_locus_all_accessions_${flanking_nucleotides}.fasta"

# follow naming scheme supported by pggb 
# [sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]

while read -r col1 col2; do
    ID=$(echo "${col1}" | sed "s/_//")
    echo "${ID}"
    sed "s/>/>${ID}#1#/" ${col2} >> $genomes
done < "${input_file}"

# index fasta
cd /path/to/pggb/
samtools faidx $genomes
 
number_of_haplotypes=$(wc -l ${genomes}.fai | awk '{print $1}')

pggb -i ${genomes} -n $number_of_haplotypes -p 95 -s 20000
