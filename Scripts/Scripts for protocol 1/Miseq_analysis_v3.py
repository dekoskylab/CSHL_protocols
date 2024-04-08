import os
import sys
import re
import subprocess
f_name=sys.argv[1]
len_gene=sys.argv[2]
template=sys.argv[3]
mm=re.search('_L001',f_name).start()
print (mm)
expt_name=f_name[0:mm]
print (expt_name)
VDJaa_fname=expt_name+'_VDJaa.txt'
print (VDJaa_fname)
cmd1='wc -l '+f_name
output=os.system(cmd1)

#Print amino acid sequences from IgBlast output files

cmd2='bash tab_to_fasta.sh ' + f_name + " " +  VDJaa_fname
print(cmd2)
os.system(cmd2)

cmd_tr = 'perl translate.pl -f 1 ' + VDJaa_fname + "_nt > " + VDJaa_fname + "_tr"
cmd_cl = "fasta_formatter -i " + VDJaa_fname + "_tr | grep -v '^>' | sed 's/.*VLA//g' | sed 's/ASTG.*//g' > " + VDJaa_fname + "_ff"

os.system(cmd_tr)
os.system(cmd_cl)

#Split files into smaller chunks
cmd3='split -l 100000 '+VDJaa_fname+"_ff"+" "+VDJaa_fname+"_ff"+"_"
os.system(cmd3)
cmd4='ls | grep '+VDJaa_fname+"_ff"+'_'
file_count=os.popen(cmd4).read()
files=[]
files=file_count.split()
print (files)

for i in range (len(files)):
        command="python jobsubmit.py "+files[i]+" "+len_gene+ " "+template
        cmd5="echo '"+command+"' >temp_command"
        print (cmd5)
        os.system(cmd5)
        cmd6="cat blank_jobsubmit temp_command"+">"+"temp_submit"
        os.system(cmd6)
        cmd7="sbatch temp_submit"
        os.system(cmd7)

cmd="python sorting.py "+f_name+" "+len_gene+ " " + template
cmd8="echo '"+cmd+"' >temp"
os.system(cmd8)
print (cmd8)
cmd9="cat blank_jobsubmit temp "+"> "+VDJaa_fname+"_ff-"+"sort_blast"
print (cmd9)
os.system(cmd9)
