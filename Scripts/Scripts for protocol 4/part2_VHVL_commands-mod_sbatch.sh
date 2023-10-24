#This file continues the VHVL_commands.sh workflow
#It pulls the necessary variables from that file, submits for parralel annotation, and also submits dependent jobs to run after that for follow-up analysis

	file1=$1
	file2=$2
	prefix=$3
	barcodes=$4

	#submit for parallel annotation of fasta files via IgBlast scripts.  Will record all the submitted job IDs in a file _queuelist.txt for later use
	rm $prefix"_queuelist.txt" 
    echo "Running blast"
	bash parallel_jobsubmit_v2-mod_sbatch.sh $prefix"_commands_part1" $prefix
	cat $prefix"_queuelist.txt" | tr -d '\n' | tr -d '\r' | sed 's/^,/:/' | sed 's/,/:/g;s/:$//' > $prefix"_queuelist2.txt"
	QL2=$(cat $prefix"_queuelist2.txt")
    echo $QL2
    echo "queuelist2 QL2" $QL2
	sleep 3
	
	#compile_igblast_isotype.sh will compile all the igblast isotype files into one
	rm $prefix"_queuelist.txt" 
    echo "Compiling blast results"
	cat blank_jobsubmit | awk -v prefix=$prefix '{print}; END {print "bash compile_igblast_isotype.sh " prefix}' > $prefix"_temp2_jobsubmit"
	sbatch --dependency=afterok:$QL2  $prefix"_temp2_jobsubmit" | awk '{print $4 ","}' >> $prefix"_queuelist.txt" 
	cat $prefix"_queuelist.txt" | tr -d '\n' | tr -d '\r' | sed 's/^,/:/' | sed 's/,/:/g;s/:$//' > $prefix"_queuelist2.txt"
	QL2=$(cat $prefix"_queuelist2.txt")
    echo $QL2
    echo "queuelist QL2" $QL2
	sleep 6
	
	#Run part2 and part 3
	rm $prefix"_queuelist.txt"
    echo "Pairing and Analysis"
	cat blank_jobsubmit | awk -v prefix=$prefix '{print}; END {print "bash CDR3motif_search_part2_v2.1.sh " prefix "_ " prefix "_igblast_isotype_all.txt\nbash CDR3motif_search_analysis_v3.2.sh " prefix "_ " prefix "_CDR3_nt_pairs_over1read.txt"}' > $prefix"_temp3_jobsubmit"
	sbatch --dependency=afterok:$QL2  $prefix"_temp3_jobsubmit" | awk '{print $4 ","}' >> $prefix"_queuelist.txt" 
	cat $prefix"_queuelist.txt" | tr -d '\n' | tr -d '\r' | sed 's/^,/:/' | sed 's/,/:/g;s/:$//' > $prefix"_queuelist2.txt"
	QL2=$(cat $prefix"_queuelist2.txt")
    echo $QL2
    echo "queuelist2 QL2" $QL2
	
	
	#Note:  If part2 and part3 above do not finish (i.e. massive dataset and clustering proceeds slowly) can modify this script to split into two different ones down the road.
	
