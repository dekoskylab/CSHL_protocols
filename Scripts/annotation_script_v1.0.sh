#!/bin/bash

#First loop (sorting loop): is to get all input files (igblast_isotype_all2.txt) sorted in this layout ( Sequence "tab" prevelence "tab" line-count)
#Input= raw data ($igblast_isotype.txt) &&&& output = sorted files ($FILENAMEsorted.txt) 
#FILELIST contains tab-separated file descriptiosn (no spaces) in Column 1 followed by "tab" then the actual filename in Column 2 followed by "tab"
# $1 is the variable for annotation_job_scripts.v.x
#Run the command as : bash annotation_script_v1.1.sh annotation_job_names.txt VL+ unsorted file

FILELIST=$1

rm sortedfiles.txt
SORTEDLIST=sortedfiles.txt
MATCH=$2

while read line
do
	#Report the file description
	DESCR=$(echo "$line" | awk -F "\t" '{print $1}')
	#Report the file itself
	FILENAME=$(echo "$line" | awk -F "\t" '{print $2}')
	
	echo "$DESCR"
	
	echo "$FILENAME"
	
	echo "PASS"
	
	COUNT=$(awk -F "\t" '{print $2}' $FILENAME | awk '{sum+=$1}; END {print sum}')	
	
    echo "COUNT=$COUNT"
    	
	awk -F "\t" '{print $0}' $FILENAME | sort -k 3,3 -u | sed 's/^ *//g' | sed "s/ /$(echo -e "\t")/g" | awk -F "\t" -v tempcount=$COUNT '{print $3 "\t" $2/tempcount "\t" $2}' > "$FILENAME"sorted.txt
	
	echo "$FILENAME"
		
	SORTEDFILE="$FILENAME"sorted.txt
	echo "$SORTEDFILE" 
	echo "$SORTEDFILE" >> sortedfiles.txt
	
	echo
    
    sleep 5

	
done < $FILELIST

#End of first loop: each file gets sorted and named as "$FILENAME"sorted.txt
	
#Second loop (joining loop): to get all input files (sorted files) joined

#SORTEDLIST ( contains all sorted files from loop 1) in column 1
# Input = ($SORTEDLIST) &&& output = ($COUNTER_jon.txt)

echo > counterfile.txt

cat $MATCH | sort -k 3,3 -u | awk -F "\t" '{print $3 "\t" $4 "\t" $5 "\t" $6}' > Seq_list.txt

JOINFILE1=Seq_list.txt

while read line
do

	#Report the file description
	DESCR=$(echo "$line" | awk -F "\t" '{print $1}')
	
	JOINFILE2=$(echo "$line" | awk -F "\t" '{print $1}')
	
	COUNTER=$(wc -l counterfile.txt | awk '{print $1}')
	echo "Loop round $COUNTER"
	echo "Join file 1 is $JOINFILE1"
	echo "Join file 2 is $JOINFILE2"
	echo

	#Print paired lines to tempfile1
	join -t "$(echo -e "\t")" $JOINFILE1 $JOINFILE2 > tempfile1.txt
	#Print unpaired lines as 0's to tempfile2
	join -v 1 -t "$(echo -e "\t")" $JOINFILE1 $JOINFILE2 | awk -F "\t" '{print $0 FS "0" FS "0"}' | sed 's/^ *//g' | sed "s/ /$(echo -e "\t")/g" > tempfile2.txt
	#Concatenate tempfile1 and tempfile2 and sort to make the final combined file for this round
	cat tempfile1.txt  tempfile2.txt | sort -k 1,1 > round"$COUNTER"_join.txt
	
	JOINFILE1=round"$COUNTER"_join.txt
	echo >> counterfile.txt
	
done < $SORTEDLIST


#Third loop: calculating enrichment ratios cross the rounds ( we replace zeros in column $2 with dashes to avoid having an error due to division by zero)

JOINFILE=round19_join.txt

echo "awk math starting"
echo

awk -F "\t" '{if($5==0){print $0}else{print $0 "\t" $7/$5 "\t" $9/$5 "\t" $11/$5 "\t" $13/$5 "\t" $15/$5 "\t" $17/$5 "\t" $19/$5 "\t" $21/$5 "\t" $23/$5 "\t" $25/$5 "\t" $27/$5 "\t" $29/$5 "\t" $31/$5 "\t" $33/$5 "\t" $35/$5 "\t" $37/$5 "\t" $39/$5 "\t" $41/$5}}' $JOINFILE > Enrich_R.txt

echo "awk math finished"
echo 

echo "Header starting"
echo 

cat annotation_job_names_NNK.txt | awk -F "\t" '{print $1}' | tr '\n' '\t' | awk -F "\t" '{print "Description" "\t"  "\t" "\t" "\t" $1 "\t\t" $2 "\t\t" $3 "\t\t" $4 "\t\t" $5 "\t\t" $6 "\t\t" $7 "\t\t" $8 "\t\t" $9 "\t\t" $10 "\t\t" $11 "\t\t" $12 "\t\t" $13 "\t\t" $14 "\t\t" $15 "\t\t" $16 "\t\t" $17 "\t\t" $18 "\t\t" $19 "\t\t" $20}' > Header.txt 

cat annotation_job_names_NNK.txt | awk -F "\t" '{print $2}' | tr '\n' '\t' | awk -F "\t" '{print "Filename" "\t" "\t" "\t" "\t" $1 "\t\t" $2 "\t\t" $3 "\t\t" $4 "\t\t" $5 "\t\t" $6 "\t\t" $7 "\t\t" $8 "\t\t" $9 "\t\t" $10 "\t\t" $11 "\t\t" $12 "\t\t" $13 "\t\t" $14 "\t\t" $15  "\t\t" $16 "\t\t" $17 "\t\t" $18 "\t\t" $19 "\t\t" $20}' >> Header.txt

cat annotation_job_names_NNK.txt | awk -F "\t" '{print $1}' | tr '\n' '\t' | awk -F "\t" '{print "Data_type" "\t" "Original" "\t" "Substituted" "\t" "Position" "\t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "Fraction \t" "Reads \t" "ER_" $2 "\t" "ER_" $3 "\t" "ER_" $4 "\t" "ER_" $5 "\t" "ER_" $6 "\t" "ER_" $7 "\t" "ER_" $8 "\t" "ER_" $9 "\t" "ER_" $10 "\t" "ER_" $11 "\t" "ER_" $12 "\t" "ER_" $13 "\t" "ER_" $14 "\t" "ER_" $15 "\t" "ER_" $16 "\t" "ER_" $17 "\t" "ER_" $18 "\t" "ER_" $19 "\t" "ER_" $20}' >> Header.txt

cat Header.txt Enrich_R.txt > BM_vFP_NNK.txt
echo "Header finished"
echo


#Rember to insert spaces at the end of the file name! othewise bash will not recognise the last file 

#Ignore round (file#+1) because it comes from joing last round with sortedfile.txt
 
