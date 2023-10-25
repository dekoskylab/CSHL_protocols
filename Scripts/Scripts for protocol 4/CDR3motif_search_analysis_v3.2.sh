#!/bin/bash
#
#Process CDR3 motif nt results file _CDR3_nt_pairs.txt via clustering
#Input:  CDR3pairing_analysis_v1.0.sh Expt_prefix EXPT_CDR3_nt_pairs.txt

EXPTNAME=$1
RAW_CDR3ntpairs=$2

rm "$EXPTNAME"unique_pairs_over1read.txt

echo "$(date +%b-%d-%H:%M:%S)   Preparing files..."

#Generating a list of HC sequences observed in pairs >=2 times
awk -v exptname="$EXPTNAME" 'BEGIN {i=1}; {if ($1>1) print ">" exptname "_completepairs_CDR-H3_Rank_" i "_\n" $2; i=i+1}' "$RAW_CDR3ntpairs" > "$EXPTNAME"heavynt_junctions_over1read.fasta



echo "$(date +%b-%d-%H:%M:%S)   96% clustering the CDR-3 junction sequences..."
usearch -cluster "$EXPTNAME"heavynt_junctions_over1read.fasta -w 4 --maxrejects 0 -usersort --iddef 2 --nofastalign -id 0.96 -minlen 11 -uc "$EXPTNAME"results.uc -seedsout "$EXPTNAME"seeds.uc

fasta_formatter -i "$EXPTNAME"seeds.uc -o "$EXPTNAME"seeds.uc_formatted
#generate list of clonotypes & frequencies based on unique sequencees
awk '/>/{getline;print}' "$EXPTNAME"seeds.uc_formatted > "$EXPTNAME"uniqueHCntseqs_over1read.txt

sort -k 1,1 "$EXPTNAME"uniqueHCntseqs_over1read.txt > "$EXPTNAME"uniqueHCntseqs_over1read.txt_sorted
sort -t $'\t' -k 2,2 -u  "$RAW_CDR3ntpairs" > "$RAW_CDR3ntpairs"_sorted 


join -t $'\t' -j1 1 -j2 2 -o 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 "$EXPTNAME"uniqueHCntseqs_over1read.txt_sorted "$RAW_CDR3ntpairs"_sorted | sort -t $'\t' -k 1,1 -nr > "$EXPTNAME"unique_pairs_over1read.txt


#Filter unique_pairs for seqs with >=2 reads
#awk '{if($1>1)print}' "$EXPTNAME"unique_pairs.txt | sort -n -r > "$EXPTNAME"unique_pairs_over1read.txt


#Shortclust workflow to account for shorter CDR-H3s
	awk 'BEGIN {temp="AAAAAAAAAAAAAAAAAAAAAAAA"}; {pad=(26-length($2)); if(length($2)<26) print ">"$2 "_pad_" pad "\n" $2 substr(temp,1,pad)}' "$EXPTNAME"unique_pairs_over1read.txt > "$EXPTNAME"temp.fasta
	usearch -cluster "$EXPTNAME"temp.fasta -w 2 --maxrejects 0 -usersort --queryalnfract 1.0 --targetalnfract 1.0 --iddef 2 --nofastalign -id 0.96 -minlen 11 -uc "$EXPTNAME"shortresults.uc -seedsout "$EXPTNAME"shortseeds.uc
	awk '{if(length($2)>=26)print}' "$EXPTNAME"unique_pairs_over1read.txt > "$EXPTNAME"longer.temp
	awk '{if(length($2)<26)print}' "$EXPTNAME"unique_pairs_over1read.txt | sort -k 2,2 > "$EXPTNAME"shorter.temp
	grep ^C "$EXPTNAME"shortresults.uc | awk '{print $9}' | sed 's/_pad_/\t/g' | awk '{print $1}' | sort > "$EXPTNAME"temp2.txt
	join -t $'\t' -j1 2 -j2 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.20 1.21 1.22 1.23 1.24 1.25 "$EXPTNAME"shorter.temp "$EXPTNAME"temp2.txt > "$EXPTNAME"temp3.txt
	cat "$EXPTNAME"temp3.txt "$EXPTNAME"longer.temp | sort -n -r | sed 's/\r/\n/g' > "$EXPTNAME"temp.txt
	date


#Convert to amino acid seqs, compile final files and FASTA files with CDR-H3 for unique pairs and top heavy/light filtered pairs

awk -v exptname="$EXPTNAME" -F "\t"  'BEGIN {i=1}; {print ">" exptname "unique_final_CDR-H3_Rank_" i "_" $1 "_reads_" $2 "_" $3 "_"  $4 "_" $5 "_" $6 "_" $7 "_" $8 "_" $9 "_" $10 "_" $11 "_" $12 "_" $13 "_" $14 " \n" $2; i=i+1}' "$EXPTNAME"temp.txt | perl translate.pl -f 1 | grep -v '>' > "$EXPTNAME"H3aatemp.txt

awk -v exptname="$EXPTNAME" -F "\t"  'BEGIN {i=1}; {print ">" exptname "unique_final_CDR-H3_Rank_" i "_" $1 "_reads_" $2 "_" $3 "_"  $4 "_" $5 "_" $6 "_" $7 "_" $8 "_" $9 "_" $10 "_" $11 "_" $12 "_" $13 "_" $14 " \n" $3; i=i+1}' "$EXPTNAME"temp.txt | perl translate.pl -f 1 | grep -v '>' > "$EXPTNAME"L3aatemp.txt

paste "$EXPTNAME"H3aatemp.txt "$EXPTNAME"L3aatemp.txt "$EXPTNAME"temp.txt | awk -F "\t" '{print $3 "\t" $4 "\t" $5 "\t" $1 "\t" $2 "\t" $6 "\t" $7 "\t" $8 "\t" (length($1)-2)  "\t" $10 "\t" $11 "\t" (length($2)-2) "\t" $13 "\t" $14}' > "$EXPTNAME"unique_final.txt

awk -v exptname="$EXPTNAME" -F "\t"  'BEGIN {i=1}; {print ">" exptname "unique_final_CDR-H3_Rank_" i "_" $1 "_reads_" $2 "_" $3 "_"  $4 "_" $5 "_" $6 "_" $7 "_" $8 "_" $9 "_" $10 "_" $11 "_" $12 "_" $13 "_" $14 " \n" $2; i=i+1}' "$EXPTNAME"unique_final.txt > "$EXPTNAME"unique_CDRH3_final.fna

awk -v exptname="$EXPTNAME" -F "\t"  'BEGIN {i=1}; {print ">" exptname "unique_final_CDR-H3_Rank_" i "_" $1 "_reads_" $2 "_" $3 "_"  $4 "_" $5 "_" $6 "_" $7 "_" $8 "_" $9 "_" $10 "_" $11 "_" $12 "_" $13 "_" $14 " \n" $4; i=i+1}' "$EXPTNAME"unique_final.txt > "$EXPTNAME"unique_CDRH3_final.faa







echo "$(date +%b-%d-%H:%M:%S)   Job complete."

echo 'Total number of paired reads:'
awk '{sum+=$1}; END {print sum}' "$RAW_CDR3ntpairs"

echo 'Recovered unique pairs:'
#wc -l "$EXPTNAME"unique_pairs.txt
wc -l "$EXPTNAME"unique_final.txt
