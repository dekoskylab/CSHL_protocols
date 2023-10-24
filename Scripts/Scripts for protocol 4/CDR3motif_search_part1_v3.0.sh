#!/bin/bash
#
#Process two FASTA files R1 and R2 for VH:VL pairs using CDR3 motif
#Input:   CDR3motif_search.sh R1_fasta_file R2_fasta_file barcodes_file
#barcodes_files :  mouse_barcodes.txt  human_barcodes.txt  rhesus_barcodes.txt

R1fasta=$1
R2fasta=$2
barcodes=$3

date; perl igblast_script.pl "$R1fasta" "$R2fasta" "$barcodes"; date

