#!/bin/bash

# This script is used to find and extract the s-locus in assemblies. Additionally, reversing of S-locus contigs is done with this script:


# Set working directories and input files
work_dir="/path/to/working/directory"
accessions_file="${work_dir}/Accessions_slocus_synteny_250000_flanks.txt"
ARK_B80="${work_dir}/S_locus/S_locus_flanks.fasta"
input_fasta="${work_dir}/All_s_locus_contigs_synteny.fasta"
output_fasta="${work_dir}/All_s_locus_contigs_corrected.fasta"
reference="${work_dir}/contig/used/as/strand/reference_s_locus_contig.fasta"
paf_output="${work_dir}/contigs_against_reference.paf"

# Function to extract S-locus contigs for each accession
extract_s_locus_contigs() {
  accession=$1

  # Define directories based on the accession
  accession_dir="${work_dir}/blast/${accession}"
  mkdir -p "$accession_dir"
  final_assembly="${work_dir}/assemblies/${accession}/intermediate_files_directories/MABS_hifiasm/mmseqs_yeast/${accession}_MABS_assembly_contamination_free.fasta"

  # Check if assembly exists
  if [ ! -f "$final_assembly" ]; then
    echo "Assembly for $accession not found, skipping." >&2
    return
  fi

  # Create BLAST database
  makeblastdb -in "$final_assembly" -dbtype nucl -out "${accession_dir}/${accession}_assembly_db"

  # Run BLAST search
  blastn -query "$ARK_B80" -db "${accession_dir}/${accession}_assembly_db" -out "${accession_dir}/${accession}_ARK_B80_match.txt" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -num_threads 4

  # Filter BLAST results for high identity matches
  awk '$3 >= 96' "${accession_dir}/${accession}_ARK_B80_match.txt" > "${accession_dir}/${accession}_ARK_B80_match_filtered.txt"
  
  # Extract contig of interest
  contig_of_interest=$(awk '{print $2}' "${accession_dir}/${accession}_ARK_B80_match_filtered.txt" | sort | uniq)
  if [ -z "$contig_of_interest" ]; then
    echo "No valid contig found for $accession."
    return
  fi

  # Extract contig sequence from assembly
  ID=$(echo "$accession" | sed "s/_//")
  grep -A 1 "$contig_of_interest" "$final_assembly" | sed "s/>/>${ID}#1#/" > "${accession_dir}/${accession}_s_locus_contig.fasta"

  # Clean up temporary files
  rm "${accession_dir}/${accession}_assembly_db"* 
}

# Export function to parallel for execution
export -f extract_s_locus_contigs
export ARK_B80

# Parallel execution for all accessions
num_jobs=$(wc -l < "$accessions_file")
cat "$accessions_file" | parallel -j "$num_jobs" extract_s_locus_contigs {}

# Concatenate all contigs into one file
> "$input_fasta"
while read -r accession; do
  cat "${work_dir}/blast/${accession}/${accession}_s_locus_contig.fasta" >> "$input_fasta"
done < "$accessions_file"

# Correct strand orientation using minimap2 and sequence reversal if necessary
> "$output_fasta"
while read -r accession; do
  accession_id=$(echo "$accession" | sed "s/_//")
  echo "Processing $accession_id"

  minimap2 -x asm10 "$reference" "$input_fasta" > "$paf_output"

  if [ "$accession_id" = "ES03014" ]; then
    echo "Skipping reference $accession_id"
    continue
  fi

  # Check strand orientation using PAF file
  sign=$(grep "$accession_id" "$paf_output" | head -n 2 | awk '{print $5}' | uniq)
  contig=$(grep -A 1 "$accession_id" "$input_fasta")

  if [ "$sign" = "-" ]; then
    echo "Contig $accession_id is on the reverse strand. Reverse complementing..."
    echo "$contig" | seqtk seq -r - >> "$output_fasta"
  else
    echo "Contig $accession_id is on the forward strand. Keeping as is..."
    echo "$contig" >> "$output_fasta"
  fi

  # Save corrected contig to individual file
  echo "$contig" > "${work_dir}/blast/${accession}/${accession}_s_locus_contig_corrected.fasta"
done < "$accessions_file"

echo "All contigs processed. Final corrected contigs saved to $output_fasta"
