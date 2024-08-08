#USAGE:  bash parallel_jobsubmit.sh commandlist prefix
#Requires template submission file blank_jobsubmit

OLD_IFS=$IFS

IFS=$'\n'
file=$1
prefix=$2

for line in $(cat $file);
	do
        echo $line
		echo $line > $prefix"_temp_command"
		IFS=$OLD_IFS
#		echo $file
		cat blank_jobsubmit $prefix"_temp_command" > $prefix"_temp_jobsubmit"
#		command=$(echo $line | awk '{print $1}')
#		file1=$(echo $line | awk '{print $2}')
#		file2=$(echo $line | awk '{print $3}')
        sbatch $prefix"_temp_jobsubmit" | awk '{print $4 ","}' >> $prefix"_queuelist.txt"
		sleep 2
	done

IFS=$OLD_IFS
