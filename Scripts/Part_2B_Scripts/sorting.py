import os
import sys
import re
import subprocess
f_name=sys.argv[1]
len_gene=sys.argv[2]
template=sys.argv[3]
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
        os.system(cmd1)

cmd1="ls "+"| grep "+expt_name+" | grep filtered_aa | xargs  cat | grep -v 'Z\|X' >> "+VDJaa_fname+"_total_aa.txt"
os.system(cmd1)

cmd2="cat "+VDJaa_fname+"_total_aa.txt"+" | sort | uniq -c >"+VDJaa_fname+"_unique_aa_sorted.txt"
os.system(cmd2)

cmd3="cat "+VDJaa_fname+"_unique_aa_sorted.txt | sort -nr -k1 >"+VDJaa_fname+"_unique_aa_sorted_reversed.txt"
os.system(cmd3)

#making database
cmd4="cat "+VDJaa_fname+"_unique_aa_sorted_reversed.txt >>"+VDJaa_fname+"_final_aa_sorted.txt"
os.system(cmd4)

cmd5="cat "+VDJaa_fname+"_final_aa_sorted.txt | awk '"+"BEGIN {i=1};"+"{print "+'">Seq"i'+'"_"$1'+";i=i+1;print $2}'"+">"+VDJaa_fname+"_aa_database.fasta"
os.system(cmd5)

#usearch

cmd3="usearch -query "+ template +" -db "+ VDJaa_fname+"_aa_database.fasta " +"-queryalnfract 1 -targetalnfract 1 --id 0.90 -userfields query+target+id+trow+gaps  --userout "+VDJaa_fname+"_user_out-AA --global  -nousort -blastout "+VDJaa_fname+"_aln_out-AA"
print(cmd3)
os.system(cmd3)

cmd3_1="head -1 "+VDJaa_fname+"_user_out-AA > user_out"
os.system(cmd3_1)

cmd4="cat "+VDJaa_fname+"_user_out-AA |  awk"+" '{if ( $5 == 0 )"+"print $1"+'"\\t"'+"$2"+'"\\t"'+"$3"+'"\\t"'+"$4}'"+" >> user_out"
print (cmd4)
os.system(cmd4)

cmd5="mv "+VDJaa_fname+"_user_out-AA "+VDJaa_fname+"_user_out-AA_old"
os.system(cmd5)

cmd6="mv user_out "+VDJaa_fname+"_user_out-AA"
os.system(cmd6)
