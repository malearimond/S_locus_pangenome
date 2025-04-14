#!/bin/bash

assembly="/path/to/input.fasta"
outname="Output_name"
outdir="/path/to/outdir"

earlGrey -g ${assembly} -s ${outname} -o ${outdir} -t 32
