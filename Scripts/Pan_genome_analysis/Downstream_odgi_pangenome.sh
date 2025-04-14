#!/bin/bash

# This script is used to perform downstream analysis of generated pan-genome graph
graph="/path/to/pggb/graph.gfa"

## paths ##
 # Which paths with which coordinates do we see in the graph
odgi paths -i ${graph} -L -l           

#### Visualization ####
## 1D including inserted genes ##
bed_file_of_genes="ARK_B80_ES03_014.bed"
odgi procbed -i ${graph} -b ${bed_file_of_genes} > ${bed_file_of_genes}.adj.bed
odgi inject -i ${graph} -b ${bed_file_of_genes}.adj.bed -o ${bed_file_of_genes}.slocus.genes.og
odgi viz -i ${bed_file_of_genes}.slocus.genes.og -o ${bed_file_of_genes}.slocus.genes.png -c 12 -w 100 -y 50 -m -B Spectral:4

## 2D: extracting path IDs and adjusting layout
# generate index:
odgi pathindex -i ${graph} -o Pathindex.xg

# define layout
odgi layout -i ${graph} -x 100 -G 100 -X Pathindex.xg -o ${layout_optimized}    # for very compact visualization
odgi layout -i ${graph} -G 9 -x 20 -X Pathindex.xg -o ${layout_optimized}     # for more flexible visualization

# without colour #
odgi draw -i ${graph} -c ${layout_optimized} -w 8000 -H 3000 -R 2 -p S_locus_${group}_nc_accessions.png

# including colour codes per paths #
# awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' tmp_colour.bed > ${bed_colours}
odgi pathindex -i ${graph} -o Pathindex.xg
odgi layout -i ${graph} -x 100 -G 100 -X Pathindex.xg -o ${layout_optimized}    # for very compact visualization
odgi layout -i ${graph} -G 9 -x 20 -X Pathindex.xg -o low_annealing_step.lay    # for more flexible visualization
odgi draw -i ${graph} -c ${layout} -w 8000 -H 3000 -R 2 -p ${output}.png -b ${bed_colours}
