#!/bin/bash

repeat_lib="/path/to/lib.fa.strained"
dir="/path/tou/output/masked/"
fasta="path/to/.fasta"
RepeatMasker --pa 1 -a -xsmall -gccalc -excln -gff -s -lib $repeat_lib -dir $dir $fasta
