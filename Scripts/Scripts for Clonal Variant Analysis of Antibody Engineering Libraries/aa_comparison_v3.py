import re
import sys
import os
import difflib
import string
import numpy as np
f_name=sys.argv[1]
len_gene=sys.argv[2]
mm=re.search('_igblast',f_name).start()
expt_name=f_name[0:mm]
result_fname=expt_name+'_VDJaa.txt'+'_user_out-AA'
FH=open(result_fname,'r')
data=FH.readlines()

#Wt_query_HC

for i in range (0,len(data)):
        mm=re.match('(^.*\tSeq\d+_\d+\t100.0)(\s+.+)',data[i])
        if mm:
                inf=data[i].split("\t")
                subject_seq=inf[3].strip()
                wt_details=inf[1].split('_')
                wt_count=wt_details[1]
                wt_seq_detail=i
                exit

fname1=expt_name+'_result_aa-comp.txt'

header = "Sequence name" +  "\t" +  "Number of reads" + "\t" +  "Amino Acid sequence" + "\t" + "Template residue" + "\t" + "Mutated residue" + "\t" + "Position"

with open(fname1, 'w') as f:
    f.write(header + "\n")
    f.close

FW=open(fname1,'a')


A_count=0;C_count=0;D_count=0;E_count=0;F_count=0;G_count=0;H_count=0;I_count=0;K_count=0;L_count=0;M_count=0;N_count=0;P_count=0;Q_count=0;R_count=0
S_count=0;T_count=0;V_count=0;W_count=0;Y_count=0; X_count=0
gene_map=[]

cmd1='cat '+expt_name+'_VDJaa.txt'+'_total_aa.txt'+" | wc -l"
total_reads=os.popen(cmd1).read()
total_reads=total_reads.rstrip("\n")


d={}
for x in range (1,int(len_gene)+1):
        d["aa{0}".format(x)]=0

row_list=[]
for m in range (0,int(len_gene)+1):
        name=int(m)+1
        row_list.append(name)

aa_mutation_count=[]

for i in range (1,len(data)):
        if (i!=wt_seq_detail):
                inf1=data[i].split("\t")
                target_id=inf1[1]
                query_seq=inf1[3].strip()
                for j in range(0,len(subject_seq)):
                        if query_seq[j]!=subject_seq[j]:
                                no_reads=target_id.split('_')
                                reads=no_reads[1]
                                reads=int(reads)
                                string=target_id+"\t"+str(reads)+"\t"+query_seq+"\t"+subject_seq[j]+"\t"+query_seq[j]+"\t"+str(int(j+1))+"\n"
                                for i in range (1,int(len_gene)+1):
                                        if (int(j+1)==i):
                                                name='aa'+str(int(j+1))
                                                d[name] += int(reads)

                                for k in range (1,int(len_gene)+1):
                                        if (int(j+1)==k):
                                                name='aa'+str(int(j+1))
                                                to_add=name+'_'+str(int(reads))+"_"+query_seq[j]
                                                aa_mutation_count.append(to_add)

                                gene_map.append(int(j+1))
                                FW.writelines (string)
                                if query_seq[j]=='A':
                                        A_count +=reads
                                elif query_seq[j]=='C':
                                        C_count+=reads
                                elif query_seq[j]=='D':
                                        D_count+=reads
                                elif query_seq[j]=='E':
                                        E_count+=reads
                                elif query_seq[j]=='F':
                                        F_count+=reads
                                elif query_seq[j]=='G':
                                        G_count+=reads
                                elif query_seq[j]=='H':
                                        H_count+=reads
                                elif query_seq[j]=='I':
                                        I_count+=reads
                                elif query_seq[j]=='K':
                                        K_count+=reads
                                elif query_seq[j]=='L':
                                        L_count+=reads
                                elif query_seq[j]=='M':
                                        M_count+=reads
                                elif query_seq[j]=='N':
                                        N_count+=reads
                                elif query_seq[j]=='P':
                                        P_count+=reads
                                elif query_seq[j]=='Q':
                                        Q_count+=reads
                                elif query_seq[j]=='R':
                                        R_count+=reads
                                elif query_seq[j]=='S':
                                        S_count+=reads
                                elif query_seq[j]=='T':
                                        T_count+=reads
                                elif query_seq[j]=='V':
                                        V_count +=reads
                                elif query_seq[j]=='W':
                                        W_count+=reads
                                elif query_seq[j]=='Y':
                                        Y_count+=reads
                                else:
                                        X_count+=reads


FW.close()
