#!/bin/bash
#SBATCH --partition=sixhour
#SBATCH --ntasks-per-node=2
#SBATCH --nodes=2
#SBATCH --constraint=ib
#SBATCH --mem-per-cpu=40gb
#SBATCH --time=0-06:00:00
#SBATCH --output=mpi_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ahmed.fahad@ku.edu

while read line
do
        f=$(echo "$line")
        file1=$(echo "$line" | awk '{print $1}')
        file2=$(echo "$line" | awk '{print $2}')
        prefix=$(echo "$line" | awk '{print $3}')
        #barcodes=$(echo "$line" | awk '{print $4}')
        barcodes=$(echo "human_barcodes.txt")
        echo $file1
        echo $file2
        echo $barcodes

        #Send command for quality filtering, calling qualityfilter_script_human_v1.0.sh.  SPECIFY THE BARCODES REQUIRED USING THE INPUT BARCODES VARIABLE AT THE TOP
        rm $prefix"_queuelist.txt"
        cat large_jobsubmit | awk -v file1=$file1 -v file2=$file2 -v prefix=$prefix -v barcodes=$barcodes '{print $0}; END {print "bash qualityfilter_script_human_v1.0.sh " file1 " " file2 " " prefix " " barcodes}' > $prefix"QFjobsubmit"
        sbatch $prefix"QFjobsubmit" | awk '{print $4 ","}' >> $prefix"_queuelist.txt"
        echo "Command ran successfully"
        cat $prefix"_queuelist.txt" | tr -d '\n' | tr -d '\r' | sed 's/^,/:/' | sed 's/,/:/g;s/:$//' > $prefix"_queuelist2.txt"
        QL2=$(cat $prefix"_queuelist2.txt")
        echo $QL2
        echo "QL2" $QL2
        sleep 3

        #Send command to begin part 2 by calling part2_VHVL_commands.sh, which will hold until the quality filtering script is complete
        echo "Making file to run part2 of pipeline"
        cat blank_jobsubmit | awk -v file1=$file1 -v file2=$file2 -v prefix=$prefix -v barcodes=$barcodes '{print $0}; END {print "bash part2_VHVL_commands-mod_sbatch.sh " file1 " " file2 " " prefix " " barcodes}' > $prefix"part2_jobsubmit"
        sbatch --dependency=afterok:$QL2 $prefix"part2_jobsubmit" | awk '{print $4 ","}' >> $prefix"_queuelist.txt"
        cat $prefix"_queuelist.txt" | tr -d '\n' | tr -d '\r' | sed 's/^,/:/' | sed 's/,/:/g;s/:$//' > $prefix"_queuelist2.txt"
        QL2=$(cat $prefix"_queuelist2.txt")
        echo $QL2

        sleep 3

done < pairlist.txt 

