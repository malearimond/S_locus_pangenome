#!/bin/bash

main_accessions=("Accession_IDs" "to" "analyze")
syri_dir="/path/to/output_dir"
chr="output_name"


# Function to rename contigs
rename_contigs() {
    local infile="$1"
    local outfile="${infile%.fasta}_renamed.fasta"
    awk '/^>/ {print ">S-locus"; next} {print}' "$infile" > "$outfile"
    echo "$outfile"
}
# Prepare genome fasta files
declare -A fasta_files
for acc in "${main_accessions[@]}"; do
      fasta_files[$acc]=$(rename_contigs "${dir}${acc}/${acc}_slocus_blast_${flanking_nucleotides}.fasta")
done


# Run nucmer comparisons
for ((i=0; i<${#main_accessions[@]}-1; i++)); do
    A="${main_accessions[$i]}"
    B="${main_accessions[$i+1]}"
    A_fa="${fasta_files[$A]}"
    B_fa="${fasta_files[$B]}"

    echo "Running nucmer for ${A} vs ${B}"
    nucmer -c 100 -b 500 -l 50 "$A_fa" "$B_fa" -p "${A}_${B}_${chr}" &
done
wait

# Run delta-filter and convert to coords
for ((i=0; i<${#main_accessions[@]}-1; i++)); do
    A="${main_accessions[$i]}"
    B="${main_accessions[$i+1]}"

    delta-filter -m -i 90 -l 100 "${A}_${B}_${chr}.delta" > "${A}_${B}_${chr}.filtered.delta" &
done
wait

for ((i=0; i<${#main_accessions[@]}-1; i++)); do
    A="${main_accessions[$i]}"
    B="${main_accessions[$i+1]}"
    show-coords -THrd "${A}_${B}_${chr}.filtered.delta" > "${A}_${B}_${chr}.filtered.coords" &
done
wait
# Create genomes.txt file
genomes="${chr}_genomes.txt"
echo -e "#file\tname\ttags" > "$genomes"
for acc in "${main_accessions[@]}"; do
    echo -e "${fasta_files[$acc]}\t${acc}\tlw:1.5" >> "$genomes"
done

# Run syri comparisons
for ((i=0; i<${#main_accessions[@]}-1; i++)); do
    A="${main_accessions[$i]}"
    B="${main_accessions[$i+1]}"

    syri -c "${A}_${B}_${chr}.filtered.coords" -d "${A}_${B}_${chr}.filtered.delta" --prefix "${A}_${B}_${chr}" -r "${fasta_files[$A]}" -q "${fasta_files[$B]}" &
done
wait

# Generate synteny plot
plot_cmd="plotsr"
for ((i=0; i<${#main_accessions[@]}-1; i++)); do
    A="${main_accessions[$i]}"
    B="${main_accessions[$i+1]}"
    plot_cmd+=" --sr \"${A}_${B}_${chr}syri.out\""
done

plot_cmd+=" --genomes \"$genomes\" -o \"Synteny_of_${chr}_final_with_genes.png\" -H 22 -W 28 -f 14 --markers markers.bed"
eval "$plot_cmd"

mkdir ${syri_dir}
mv *${chr}* $syri_dir
echo "Analysis complete!"
