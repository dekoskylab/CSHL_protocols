import os
import sys
import re
import subprocess
f_name=sys.argv[1]
len_gene=sys.argv[2]
mm=re.search('_L001',f_name).start()
print (mm)
expt_name=f_name[0:mm]
print (expt_name)
VDJaa_fname=expt_name+'_VDJaa.txt'
print (VDJaa_fname)
cmd1='wc -l '+f_name 
output=os.system(cmd1)



cmd2='cat '+f_name+ "|awk '"+'{if(length($3)>5) print $10 }'+"' >"+VDJaa_fname
#print (cmd2)
os.system(cmd2)
cmd2_1='cat '+f_name+ "|awk '"+'{if(length($4)>5) print $12 }'+"' >"+VDJaa_fname+"_addn"
os.system(cmd2_1)
cmd2_2='cat '+ VDJaa_fname+" "+ VDJaa_fname+"_addn "+" > "+VDJaa_fname+"_ff"
print (cmd2_2)
os.system(cmd2_2)


#cmd2_1='cat '+f_name+ "|awk '"+'{if(length($3)>5) print $9 }'+"' >"+VDJaa_fname+"-NT"
#os.system(cmd2_1)


cmd3='split -l 100000 '+VDJaa_fname+"_ff"+" "+VDJaa_fname+"_ff"+"_"
os.system(cmd3)
#cmd3_1='split -l 120000 '+VDJaa_fname+"-NT"+"  "+VDJaa_fname+"-NT"+"_"
#os.system(cmd3_1)

cmd4='ls | grep '+VDJaa_fname+"_ff"+'_'
file_count=os.popen(cmd4).read()
files=[]
files=file_count.split()
print (files)


for i in range (len(files)):
	command="python jobsubmit.py "+files[i]+" "+len_gene
	cmd5="echo '"+command+"' >temp_command"
	print (cmd5)
	os.system(cmd5)
	cmd6="cat blank_jobsubmit temp_command"+">"+"temp_submit"
	os.system(cmd6)
	cmd7="sbatch temp_submit"
	os.system(cmd7)

cmd="python sorting.py "+f_name+" "+len_gene
cmd8="echo '"+cmd+"' >temp"
os.system(cmd8)
print (cmd8)
cmd9="cat blank_jobsubmit temp "+"> "+VDJaa_fname+"_ff-"+"sort_blast"
print (cmd9)
os.system(cmd9)

	
	

