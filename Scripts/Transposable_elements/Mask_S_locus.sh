#!/bin/bash

repeat_lib="/path/to/lib.fa.strained"
dir="/path/tou/output/masked/"
fasta="path/to/.fasta"
RepeatMasker --pa 1 -a -xsmall -gccalc -excln -gff -s -lib $repeat_lib -dir $dir $fasta


# Process TE analysis:
output_file="repeat_summary.tsv"
echo -e "Genome_ID\tTotal_repeats\tSINE (cov [%])\tSINE (elements [#])\tLINE (cov [%])\tLINE (elements [#])\tLTR (cov [%])\tLTR (elements [#])\tDNA TE (cov [%])\tDNA TE (elements [#])\tUnclassified (cov [%])\tUnclassified (elements [#])" > $output_file

# Loop through all .tbl files
while read -r col1; do
  file="/netscratch/dep_coupland/grp_fulgione/male/danijel_s_locus/flank_0/${col1}/masked/${col1}_slocus_blast_0.fasta.tbl"
    genome_id=$col1  # Extract genome ID from filename
    total_repeats=$(grep "Total interspersed repeats:" "$file" | awk '{print $(NF-1)}')
    sine_cov=$(grep -A1 "SINEs:" "$file" | awk 'NR==1 {print $(NF-1)}')
    sine_elements=$(grep -A1 "SINEs:" "$file" | awk 'NR==1 {print $2}')
    line_cov=$(grep -A1 "LINEs:" "$file" | awk 'NR==1 {print $(NF-1)}')
    line_elements=$(grep -A1 "LINEs:" "$file" | awk 'NR==1 {print $2}')
    ltr_cov=$(grep -A1 "LTR elements:" "$file" | awk 'NR==1 {print $(NF-1)}')
    ltr_elements=$(grep "LTR elements:" "$file" | awk '{print $3}')
    dna_cov=$(grep -A1 "DNA elements:" "$file" | awk 'NR==1 {print $(NF-1)}')
    dna_elements=$(grep  "DNA elements:" "$file" | awk '{print $3}')
    unclassified_cov=$(grep -A1 "Unclassified:" "$file" | awk 'NR==1 {print $(NF-1)}')
    unclassified_elements=$(grep -A1 "Unclassified:" "$file" | awk 'NR==1 {print $2}')

    # Append results to the output file
    echo -e "$genome_id\t$total_repeats\t$sine_cov\t$sine_elements\t$line_cov\t$line_elements\t$ltr_cov\t$ltr_elements\t$dna_cov\t$dna_elements\t$unclassified_cov\t$unclassified_elements" >> $output_file
  done < /file/with/accessions.txt
