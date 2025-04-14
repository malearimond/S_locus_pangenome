#!/bin/bash

# ------------------------------
# Set working directory & inputs
# ------------------------------
work_dir="/path/to/working_dir/"
accessions="${work_dir}All_accessions_plus_OG_acc.txt"
flanking_nucleotides=0

# Define reference genes to BLAST
ARK_B80="/path/to/S_locus_flanks.fasta"

# ---------------------------------
# Function to process each accession
# ---------------------------------
process_accession () {
  accession=$1
  echo "[INFO] Processing accession: $accession"

  blast_dir="${work_dir}/blast/${accession}"
  mkdir -p "$blast_dir"
  cd "$blast_dir" || exit

  # Path to the corrected contig from script 1
  contig_fasta="${blast_dir}/${accession}_s_locus_contig_corrected.fasta"

  # Skip if the input fasta doesn't exist
  if [[ ! -f "$contig_fasta" ]]; then
    echo "[WARN] Contig file missing: $contig_fasta" >&2
    return
  fi

  # Index contig FASTA
  samtools faidx "$contig_fasta"

  # Build BLAST DB
  makeblastdb -in "$contig_fasta" -dbtype nucl -out "${accession}_assembly_db"

  # Run BLAST
  blastn -query "$ARK_B80" -db "${accession}_assembly_db" \
    -out "${accession}_ARK_B80_match.txt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -evalue 1e-5 -num_threads 4

  # Filter matches with length ≥ 300bp
  awk '$4 >= 300' "${accession}_ARK_B80_match.txt" > "${accession}_ARK_B80_match_filtered.txt"

  # Extract coordinates with correct orientation
  coords_out="${accession}_ARK_B80_coordinates.txt"
  grep 'tg' "${accession}_ARK_B80_match_filtered.txt" | awk -v out="$coords_out" '{
    if ($9 > $10) {
      print $2 "\t" $10 "\t" $9 > out
    } else {
      print $2 "\t" $9 "\t" $10 > out
    }
  }'

  # ---------------------------
  # Extract flanking FASTA
  # ---------------------------
  flank_dir="${work_dir}/flank_${flanking_nucleotides}/${accession}"
  mkdir -p "$flank_dir"
  cd "$flank_dir" || exit

  bed_file="${accession}_slocus_blast_${flanking_nucleotides}.bed"
  fasta_file="${accession}_slocus_blast_${flanking_nucleotides}.fasta"

  max_contig_value=$(awk '{print $2}' "${contig_fasta}.fai")

  awk -v flanks="$flanking_nucleotides" -v max="$max_contig_value" '{
    if (NR == 1) {
      min = $2; max_val = $3; contig = $1
    }
    if ($2 < min) min = $2;
    if ($3 > max_val) max_val = $3;
  }
  END {
    start = min - flanks;
    if (start < 0) start = 0;
    end = max_val + flanks;
    if (end > max) end = max;
    print contig "\t" start "\t" end;
  }' "$coords_out" > "$bed_file"

  # Extract FASTA using bedtools
  bedtools getfasta -fi "$contig_fasta" -bed "$bed_file" -fo "$fasta_file"

  # Log fasta location
  echo "${accession} ${fasta_file}" >> "${work_dir}/flank_${flanking_nucleotides}/OG_accessions_slocus_${flanking_nucleotides}_flanks.txt"
}

# ------------------------------------
# Export and run in parallel
# ------------------------------------
export -f process_accession
export ARK_B80
export flanking_nucleotides
export work_dir

parallel -j "$(wc -l < "$accessions")" process_accession :::: "$accessions"
