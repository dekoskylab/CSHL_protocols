import re
import sys
import os
import difflib
import string
import numpy as np
f_name=sys.argv[1]
len_gene=sys.argv[2]
mm=re.search('_L001',f_name).start()
expt_name=f_name[0:mm]
result_fname=expt_name+'_VDJaa.txt'+'_user_out-AA'
print (result_fname)
FH=open(result_fname,'r')
data=FH.readlines()

#Wt_query_HC
#Change the header based on query_aa.fasta file in re.match()

for i in range (0,len(data)):
        mm=re.match('(^.*\tSeq\d+_\d+\t100.0)(\s+.+)',data[i])
        if mm:
                print (i,data[i])
                inf=data[i].split("\t")
                subject_seq=inf[3].strip()
                wt_details=inf[1].split('_')
                wt_count=wt_details[1]
                print (wt_count,"\t",subject_seq)
                wt_seq_detail=i
                exit

fname1=expt_name+'_result_aa-comp.txt'
FW=open(fname1,'w')


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
        #print (name)

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
                        #print (target_id,str(reads),query_seq,subject_seq[j],query_seq[j],j+1)
                                string=target_id+"\t"+str(reads)+"\t"+query_seq+"\t"+subject_seq[j]+"\t"+query_seq[j]+"\t"+str(int(j+1))+"\n"
                                for i in range (1,int(len_gene)+1):
                                        if (int(j+1)==i):
                                        #print (i," ",str(reads))
                                                name='aa'+str(int(j+1))
                                                d[name] += int(reads)

                                for k in range (1,int(len_gene)+1):
                                        if (int(j+1)==k):
                                                #print (k, " ", str(reads))
                                                name='aa'+str(int(j+1))
                                                to_add=name+'_'+str(int(reads))+"_"+query_seq[j]
                                        #print (to_add)
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

print ("A = " ,A_count);#FW.write("A\t%s\n" %A_count)
print ("C = " ,C_count);#FW.write("C\t%s\n" %C_count)
print ("D = " ,D_count);#FW.write("D\t%s\n" %D_count)
print ("E = " ,E_count);#FW.write("E\t%s\n" %E_count)
print ("F = " ,F_count);#FW.write("F\t%s\n" %F_count)
print ("G = " ,G_count);#FW.write("G\t%s\n" %G_count)
print ("H = " ,H_count);#FW.write("H\t%s\n" %H_count)
print ("I = " ,I_count);#FW.write("I\t%s\n" %I_count)
print ("K = " ,K_count);#FW.write("K\t%s\n" %K_count)
print ("L = " ,L_count);#FW.write("L\t%s\n" %L_count)
print ("M = " ,M_count);#FW.write("M\t%s\n" %M_count)
print ("N = " ,N_count);#FW.write("N\t%s\n" %N_count)
print ("P = " ,P_count);#FW.write("P\t%s\n" %P_count)
print ("Q = " ,Q_count);#FW.write("Q\t%s\n" %Q_count)
print ("R = " ,R_count);#FW.write("R\t%s\n" %R_count)
print ("S = " ,S_count);#FW.write("S\t%s\n" %S_count)
print ("T = " ,T_count);#FW.write("T\t%s\n" %T_count)
print ("V = " ,V_count);#FW.write("V\t%s\n" %V_count)
print ("W = " ,W_count);#FW.write("W\t%s\n" %W_count)
print ("Y = ", Y_count);#FW.write("Y\t%s\n" %Y_count)
print ("X = ", X_count);#FW.write("X\t%s\n" %X_count)
total=A_count+C_count+D_count+E_count+F_count+G_count+H_count+I_count+K_count+L_count+M_count+N_count+P_count+Q_count+R_count+S_count+T_count+V_count+W_count+Y_count
print ("TOTAL READS =",total);#FW.write("TOTAL READS = %s"%total)
print (gene_map)

for k in range(0,int(len_gene)):
        print (k+1,"\t",gene_map.count(k+1))
FW.close()
for m, v in d.items():
        print (m,v)

#print (aa_mutation_count)

list=[[0]*20 for n in range (int(len_gene))]
for i in range (0,len(aa_mutation_count)):
        inf2=aa_mutation_count[i].split('_')
        aa_no = inf2[0]
        aa_no = aa_no[2:len(aa_no)]
        aa_reads=inf2[1]
        aa_type=inf2[2]
        #print (type(aa_no))
        for j in range (1,int(len_gene)+1):
                if (int(aa_no) == j):
                        if (aa_type == 'A'):
                                list[int(aa_no)-1][0]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'C'):
                                list[int(aa_no)-1][1]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'D'):
                                list[int(aa_no)-1][2]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'E'):
                                list[int(aa_no)-1][3]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'F'):
                                list[int(aa_no)-1][4]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'G'):
                                list[int(aa_no)-1][5]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'H'):
                                list[int(aa_no)-1][6]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'I'):
                                list[int(aa_no)-1][7]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'K'):
                                list[int(aa_no)-1][8]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'L'):
                                list[int(aa_no)-1][9]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'M'):
                                list[int(aa_no)-1][10]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'N'):
                                list[int(aa_no)-1][11]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'P'):
                                list[int(aa_no)-1][12]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'Q'):
                                list[int(aa_no)-1][13]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'R'):
                                list[int(aa_no)-1][14]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'S'):
                                list[int(aa_no)-1][15]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'T'):
                                list[int(aa_no)-1][16]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'V'):
                                list[int(aa_no)-1][17]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'W'):
                                list[int(aa_no)-1][18]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
                        if (aa_type == 'Y'):
                                list[int(aa_no)-1][19]=int(aa_reads)
                        #       print (aa_no,"\t",aa_type,"\t",aa_reads)
