#!/bin/bash

# Define paths and variables
outdir="Output_dir/of/mmseq/"
fa="/Path_to/Assembly.fasta"
db_path="/Path_to/mmseqs2/"
db="${db_path}/swissprot"

# Prepare the output directory
rm -rf "${outdir}"
mkdir -p "${outdir}" && cd "${outdir}"

# Create the MMseqs database
mmseqs createdb "${fa}" "${col1}_MABS_db"
query="${outdir}/${col1}_MABS_db"

# Perform taxonomy analysis
mmseqs taxonomy "${query}" "${db}" "${col1}_firstResult" tmp --remove-tmp-files 0
mmseqs createtsv "${query}" "${col1}_firstResult" "${col1}_taxonomyResult.tsv"

# Generate contamination report
mmseqs taxonomyreport "${db}" "${col1}_firstResult" "${col1}_contamination_report.html" --report-mode 1

# Run easy taxonomy classification
mmseqs easy-taxonomy "${fa}" "${db}" "${col1}_mmseqResult" tmp --remove-tmp-files 0

# Filter contigs for Viridiplantea (plants) and yeast
mmseqs filtertaxdb "${db}" ./tmp/latest/result "${col1}_mmseqResult_yeast" --taxon-list 4893
mmseqs filtertaxdb "${db}" ./tmp/latest/result "${col1}_mmseqResult_plants" --taxon-list 33090

# Create TSV files for scaffolds
mmseqs createtsv ./tmp/latest/query "${col1}_mmseqResult_plants" "${col1}_plants_scaffolds.tsv"
mmseqs createtsv ./tmp/latest/query "${col1}_mmseqResult_yeast" "${col1}_yeast_scaffolds.tsv"

# Generate a contamination-free FASTA
seqkit grep -f <(cut -f1 "${col1}_plants_scaffolds.tsv") "${fa}" > "$(basename "${fa%.*}")_contamination_free.fa"
seqkit grep -f <(cut -f1 "${col1}_yeast_scaffolds.tsv") "${fa}" >> "$(basename "${fa%.*}")_contamination_free.fa"

# Format the final contamination-free assembly
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' "$(basename "${fa%.*}")_contamination_free.fa" | \
tr "\t" "\n" > "${col1}_MABS_assembly_contamination_free.fasta"

# Cleanup intermediate files
rm -rf ./tmp ./*mmseqResult* ./*firstResult*
rm "$(basename "${fa%.*}")_contamination_free.fa"
