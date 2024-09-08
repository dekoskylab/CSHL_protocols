file1=$1
	file2=$2
	prefix=$3
	barcodes=$4

	echo $file1
	echo $file2
	echo $prefix
	echo $barcodes
	echo "DONE ECHO"
	
	fastq_quality_filter -Q33 -q 20 -p 50 -i $file1 -o $file1"_q20p50.fastq"
	fastq_quality_filter -Q33 -q 20 -p 50 -i $file2 -o $file2"_q20p50.fastq"
	date
	echo "$f quality filtered"
	fastq_to_fasta -Q33 -n -i $file1"_q20p50.fastq" -o $file1"_q20p50.fasta"
	fastq_to_fasta -Q33 -n -i $file2"_q20p50.fastq" -o $file2"_q20p50.fasta"
	echo "$f converted to fasta"
	awk '/>/{temp=$0; getline; print temp "\t" $0}' $file1"_q20p50.fasta" | sed 's/\/1/ 1:N:0:8/' | sort > $file1"_q20p50.fasta_temp"
	awk '/>/{temp=$0; getline; print temp "\t" $0}' $file2"_q20p50.fasta" | sed 's/\/2/ 2:N:0:8/' | sort > $file2"_q20p50.fasta_temp"
	echo "$f converted to temporary file for joining"
	join -j 1 -o 1.1 1.2 1.3 $file1"_q20p50.fasta_temp" $file2"_q20p50.fasta_temp" | awk '{print $1 " " $2 "\n" $3}' > $prefix"_r1_q20p50_filtered.fasta"
	join -j 1 -o 2.1 2.2 2.3 $file1"_q20p50.fasta_temp" $file2"_q20p50.fasta_temp" | awk '{print $1 " " $2 "\n" $3}' > $prefix"_r2_q20p50_filtered.fasta"
	echo "$f joined into filtered files"
	split -l 200000 $prefix"_r1_q20p50_filtered.fasta" $prefix"_r1_q20p50_filtered.fasta_"
	split -l 200000 $prefix"_r2_q20p50_filtered.fasta" $prefix"_r2_q20p50_filtered.fasta_"
	
		
	ls | grep $prefix"_r2_q20p50_filtered.fasta_" | grep -v temp | grep -v reads | awk '{print $1 " " $1}' | sed 's/_r2_q20p50/_r1_q20p50/' | awk -v barcodes=$barcodes '{print "CDR3motif_search_part1_v2.1.sh " $0 " " barcodes}' > $prefix"_commands_part1"