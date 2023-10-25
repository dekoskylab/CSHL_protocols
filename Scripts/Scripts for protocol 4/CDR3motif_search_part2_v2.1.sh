#!/bin/bash
#
#Process concatenated CDR3 motif file into a series of VH:VL pared raw data files and fasta files
#Input:  CDR3motif_search_part2_v2.0.sh Expt_prefix paired_igblast_isotype_file

Exptname=$1
paired_igblast_isotype_file=$2

date

echo "Compiling raw nt and aa pairings..."
awk '{if($12=="R1H") print $2 "\t" $3 "\t" length($2) "\t" length($3); else print $2 "\t" $3 "\t" length($2) "\t" length($3)}' "$paired_igblast_isotype_file" | grep -v CDRH3 | sort | uniq -c | sort -n -r > "$Exptname"CDR3_seqonly_raw_aa_pairs.txt
awk '{if($12=="R1H") print $4 "\t" $5 "\t" length($4) "\t" length($5); else print $4 "\t" $5 "\t" length($4) "\t" length($5)}' "$paired_igblast_isotype_file" | grep -v CDRH3 | sort | uniq -c | sort -n -r > "$Exptname"CDR3_seqonly_raw_nt_pairs.txt
awk '{if($12=="R1H") print $2 "\t" $3 "\t" $14 "\t" $18 "\t" $16 "\t" length($2) "\t" $15 "\t" $19 "\t" $17 "\t" length($3) "\t" $20 "\t" $21; else print $2 "\t" $3 "\t" $15 "\t" $19 "\t" $17 "\t" length($2) "\t" $14 "\t" $18 "\t" $16 "\t" length($3) "\t" $20 "\t" $21}' "$paired_igblast_isotype_file" | grep -v CDRH3 | sort | uniq -c | sort -n -r > "$Exptname"CDR3_raw_aa_pairs.txt
awk '{if($12=="R1H") print $4 "\t" $5 "\t" $14 "\t" $18 "\t" $16 "\t" length($4) "\t" $15 "\t" $19 "\t" $17 "\t" length($5) "\t" $20 "\t" $21; else print $4 "\t" $5 "\t" $15 "\t" $19 "\t" $17 "\t" length($4) "\t" $14 "\t" $18 "\t" $16 "\t" length($5) "\t" $20 "\t" $21}' "$paired_igblast_isotype_file" | grep -v CDRH3 | sort | uniq -c | sort -n -r > "$Exptname"CDR3_raw_nt_pairs.txt

date
echo "Compiling filtered list of nt and aa pairings..."
awk '{if(length($2)>1) if (length($3)>1) if ($2!=$3) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $11 "\t" $12 "\t" $13}' "$Exptname"CDR3_raw_aa_pairs.txt > "$Exptname"CDR3_aa_pairs.txt
awk '{if(length($2)>1) if (length($3)>1) if ($2!=$3) if ($4!~/IG[KL]/) if ($8!~/IGH/) if ($4!="-" && $8!="-" )  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $11 "\t" $12 "\t" $13}' "$Exptname"CDR3_raw_nt_pairs.txt > "$Exptname"CDR3_nt_pairs.txt


date
echo "Generating supplementary files..."
awk '{if($1>1)print}' "$Exptname"CDR3_aa_pairs.txt > "$Exptname"CDR3_aa_pairs_over1read.txt
awk '{if($1>1)print}' "$Exptname"CDR3_nt_pairs.txt > "$Exptname"CDR3_nt_pairs_over1read.txt


awk '{print ">CDRH3_" $1 "\n" $2}' "$Exptname"CDR3_nt_pairs_over1read.txt  | perl translate.pl -f 1 | grep -v '>' > "$Exptname"CDRH3_pairs_over1read_translated.txt

awk '{print ">CDRL3_" $1 "\n" $3}' "$Exptname"CDR3_nt_pairs_over1read.txt  | perl translate.pl -f 1  | grep -v '>' > "$Exptname"CDRL3_pairs_over1read_translated.txt

#Preparing CDR3_nt+aa_pairs_over1read.txt


paste "$Exptname"CDRH3_pairs_over1read_translated.txt "$Exptname"CDRL3_pairs_over1read_translated.txt "$Exptname"CDR3_nt_pairs_over1read.txt | awk '{print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14}' > "$Exptname"CDR3_nt+aa_pairs_over1read.txt

awk 'BEGIN {i=1}; {print ">CDRH3_rank_" i "_CDRL3pair_" $3 "_" $1 "_reads\n" $2; i=i+1}' "$Exptname"CDR3_aa_pairs_over1read.txt > "$Exptname"CDRH3_over1read.faa
awk 'BEGIN {i=1}; {print ">CDRH3_rank_" i "_CDRL3pair_" $3 "_" $1 "_reads\n" $2; i=i+1}' "$Exptname"CDR3_nt_pairs_over1read.txt > "$Exptname"CDRH3_over1read.fna

awk 'BEGIN {i=1}; {print ">CDRH3_rank_" i "_CDRL3pair_" $3 "_" $1 "_reads\n" $2; i=i+1}' "$Exptname"CDR3_nt_pairs.txt > "$Exptname"CDRH3_allreads.fna


echo "Job complete."
date
