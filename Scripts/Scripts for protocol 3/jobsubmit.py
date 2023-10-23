import os
import sys
VDJaa_fname=sys.argv[1]
VDJnt_fname=VDJaa_fname[0:len(VDJaa_fname)-3]
prefix=VDJaa_fname[-2:]
print(VDJnt_fname)
len_gene=sys.argv[2]
seq_file = sys.argv[3]

template=open(seq_file)
aa = template.readlines()

seq_temp = aa[1].rstrip('\n')

cmd3 = "awk '{print match($1,/[ACDEFGHIKLMNPQRSTVWY]" + seq_temp[1:3] + "|" + seq_temp[0] + "[ACDEFGHIKLMNPQRSTVWY]" + seq_temp[2:3] + "|" + seq_temp[0:1] + "[ACDEFGHIKLMNPQRSTVWY]" + seq_temp[3] + "|" + seq_temp[0:3] + "[ACDEFGHIKLMNPQRSTVWY]" + "/)}' " + VDJaa_fname
result=os.popen(cmd3).read()
list=[];list=result.split()
print("%s, Sequences with gene %s"%(VDJaa_fname,len(list)))
count_0=0;count_10=0;count_1_9=0;count_remain=0
for i in range (len(list)):
                if (list[i] == '10'):
                                count_10+=1; val=int(list[i])+int(len_gene)-1
                                #gene_start=(int(list[i])-1)*3+1;gene_end=val*3
                                cmd4="awk 'NR=="+str(int(i)+1)+"{print; exit}' "+VDJaa_fname+" | cut -c "+list[i]+"-"+str(val)+">>"+VDJaa_fname+"-10.txt"
                                os.system(cmd4)


                elif (list[i]== '0'):
                                count_0+=1
                elif (list[i] >= '1' and list[i] <= '9'):
                                count_1_9+=1;val=int(list[i])+int(len_gene)-1
                                #gene_start=(int(list[i])-1)*3+1;gene_end=val*3
                                cmd4="awk 'NR=="+str(int(i)+1)+"{print; exit}' "+VDJaa_fname+" | cut -c "+list[i]+"-"+str(val)+">>"+VDJaa_fname+"-1-9.txt"
                                os.system(cmd4)


                else:
                                count_remain+=1;val=int(list[i])+int(len_gene)-1
                                #gene_start=(int(list[i])-1)*3+1;gene_end=val*3
                                cmd4="awk 'NR=="+str(int(i)+1)+"{print; exit}' "+VDJaa_fname+" | cut -c "+list[i]+"-"+str(val)+">>"+VDJaa_fname+"-11.txt"
                                os.system(cmd4)



print (count_0);print(count_1_9);print (count_10);print (count_remain)
print (count_0+count_1_9+count_10+count_remain)
