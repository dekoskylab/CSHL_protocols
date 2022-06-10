import os
import sys
import re
import subprocess
f_name=sys.argv[1]
len_gene=sys.argv[2]
mm=re.search('_L001',f_name).start()
expt_name=f_name[0:mm]
VDJaa_fname=expt_name+'_VDJaa.txt'
VDJnt_fname=VDJaa_fname[0:len(VDJaa_fname)-3]
prefix=VDJaa_fname[-2:]
cmd_a="rm "+VDJaa_fname+"_final_sorted.txt"
os.system(cmd_a)
cmd_b="rm "+VDJaa_fname+"_genes_database.fasta"
os.system(cmd_b)

cmd_c="ls | grep genes | xargs rm"
os.system(cmd_c)

cmd_d=" ls | grep "+ expt_name + "| grep "+" 'filtered\|total_aa\|unique_aa_sorted\|final_aa_sorted\|aa_database\|sorted_reversed\|user_out\|aln_out'"+"  | xargs rm"
os.system(cmd_d)

len_gene=sys.argv[2]
cmd1="ls | grep '"+VDJaa_fname+"_ff"+"_[a-z][a-z]$'"
file_count=os.popen(cmd1).read()
files=[]
files=file_count.split()
print (files)


#concatinating filtered reads and sorting
for i in range (len(files)):
	cmd1="ls | grep '"+files[i]+"-10.txt\|"+files[i]+"-1-9.txt\|"+files[i]+"-11.txt'"+"| xargs cat >> "+files[i]+'_filtered_aa.txt'
	#print(cmd1)	
	os.system(cmd1)
	#cmd2="ls | grep '"+files[i]+"_NT-10.txt\|"+files[i]+"_NT-1-9.txt\|"+files[i]+"_NT-11.txt'"+"| xargs cat >> "+files[i]+'_filtered_nt.txt'	
	#os.system(cmd2)


	
cmd1="ls "+"| grep "+expt_name+" | grep filtered_aa | xargs  cat | grep -v 'Z\|X' >> "+VDJaa_fname+"_total_aa.txt"
os.system(cmd1)
#cmd1_1="ls | grep filtered_nt | xargs  cat >> "+VDJaa_fname+"_total_nt.txt"
#os.system(cmd1_1)

cmd2="cat "+VDJaa_fname+"_total_aa.txt"+" | sort | uniq -c >"+VDJaa_fname+"_unique_aa_sorted.txt" 
os.system(cmd2)
#cmd2_1="cat "+VDJaa_fname+"_total_nt.txt"+" | sort | uniq -c >"+VDJaa_fname+"_unique_nt_sorted.txt"
#os.system(cmd2_1)

cmd3="cat "+VDJaa_fname+"_unique_aa_sorted.txt | sort -nr -k1 >"+VDJaa_fname+"_unique_aa_sorted_reversed.txt"
os.system(cmd3)
#cmd3_1="cat "+VDJaa_fname+"_unique_nt_sorted.txt | sort -nr -k1 >"+VDJaa_fname+"_unique_nt_sorted_reversed.txt"
#os.system(cmd3_1)



#making database
cmd4="cat "+VDJaa_fname+"_unique_aa_sorted_reversed.txt >>"+VDJaa_fname+"_final_aa_sorted.txt"
#print(cmd4)
os.system(cmd4)
#cmd4_1="cat "+VDJaa_fname+"_unique_nt_sorted_reversed.txt >>"+VDJaa_fname+"_final_nt_sorted.txt"
#os.system(cmd4_1)

cmd5="cat "+VDJaa_fname+"_final_aa_sorted.txt | awk '"+"BEGIN {i=1};"+"{print "+'">Seq"i'+'"_"$1'+";i=i+1;print $2}'"+">"+VDJaa_fname+"_aa_database.fasta"
#print(cmd5)
os.system(cmd5)
#cmd5_1="cat "+VDJaa_fname+"_final_nt_sorted.txt | awk '"+"BEGIN {i=1};"+"{print "+'">Seq"i'+'"_"$1'+";i=i+1;print $2}'"+">"+VDJaa_fname+"_nt_database.fasta"
#os.system(cmd5_1)




#blast


cmd3="usearch -query query_aa.fasta"+" -db "+ VDJaa_fname+"_aa_database.fasta " +"-queryalnfract 1 -targetalnfract 1 --id 0.90 -userfields query+target+id+trow+gaps  --userout "+VDJaa_fname+"_user_out-AA --global  -nousort -blastout "+VDJaa_fname+"_aln_out-AA"
print(cmd3)
os.system(cmd3)

cmd3_1="head -1 "+VDJaa_fname+"_user_out-AA > user_out"
os.system(cmd3_1)

cmd4="cat "+VDJaa_fname+"_user_out-AA |  awk"+" '{if ( $5 == 0 )"+"print $1"+'"\\t"'+"$2"+'"\\t"'+"$3"+'"\\t"'+"$4}'"+" >> user_out"
print (cmd4)
os.system(cmd4)

cmd5="mv "+VDJaa_fname+"_user_out-AA "+VDJaa_fname+"_user_out-AA_old"
#print (cmd5)
os.system(cmd5)

cmd6="mv user_out "+VDJaa_fname+"_user_out-AA"
#print (cmd6)
os.system(cmd6)



'''
#cmd3_1="usearch -ublast query_nt.fasta"+" -db "+ VDJaa_fname+"_nt_database.fasta " +"-query_cov 1 -target_cov 1 --id 0.99 -userfields query+target+id  --userout "+VDJaa_fname+"_user_out-NT -alnout "+VDJaa_fname+"_aln_out-NT"
	
'''
	
