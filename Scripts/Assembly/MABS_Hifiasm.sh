# Run MABS-hifiasm including ONT data
local_busco_dataset="${directory_to_busco}/brassicales_odb10"
Hifi_reads="Path_to_hifi/.fasta"
ONT_reads="Path_to_ONT/.fasta"
output_dir="Path_to/output_directory/"
threads="Number_of_threads"

mabs_v2.28 mabs-hifiasm.py --pacbio_hifi_reads ${Hifi_reads} --ultralong_nanopore_reads ${ONT_reads} --local_busco_dataset ${local_busco_dataset} --threads ${threads} --output_folder ${output_dir}
