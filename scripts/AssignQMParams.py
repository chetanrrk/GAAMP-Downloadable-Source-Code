import os,sys

def writeQMSetting(GAAMP_RUN,QM_MEM,QM_NPROC):
    gaampScriptPath = sys.argv[0].split("AssignQMParams.py")[0]
    fin = open(gaampScriptPath+'QM-para.txt','r')
    fout = open(GAAMP_RUN+'/QM-para.txt','w')
    for line in fin:
	token=line.rstrip().split()
	if len(token)>0:
	   if token[0]=="QM_MEM":
	    	line = line.replace("QM_MEM        1GB","QM_MEM        "+QM_MEM)
	   elif token[0]=="QM_NPROC":
	    	line = line.replace("QM_NPROC      8","QM_NPROC      "+QM_NPROC)		
	print>>fout,line.rstrip()
    fout.close()
    fin.close()	


if __name__=="__main__":
    GAAMP_RUN = sys.argv[1] ## GAAMP Working Dir
    QM_MEM = sys.argv[2] ## QM Memory
    QM_NPROC = sys.argv[3] ## QM Num. Processors 	 
    writeQMSetting(GAAMP_RUN,QM_MEM,QM_NPROC)

