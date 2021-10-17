import sys, os, math
from scipy import stats
import numpy as np
dir,sampleInfo,nwise=sys.argv[1],sys.argv[2],int(sys.argv[3])
def grabInfo(sampleInfo):
 file=open(sampleInfo,"r")
 refd={}
 for ln in file:
  ln=ln.strip("\r\n").split(",")
  if ln[-1] != "": 
   if ln[-1] not in refd:
    refd[ln[-1]]=[ln[-2]]
   else:
    refd[ln[-1]].append(ln[-2])
 return refd
refd=grabInfo(sampleInfo)

filelist=[fn for fn in os.listdir(dir) if fn.startswith("BC")]
def CPMd(filelist,nwise):
 cpmd={}
 for fn in filelist:
  #file=open("/"+dir+"/"+fn,"r")
  file=open(dir+"/"+fn,"r")
  samp=fn.split("_")[1].split(".")[0]
  for y in range (8):
   file.readline()
  for ln in file:
   ln=ln.strip("\r\n").split(",")
   if ln[nwise] not in cpmd:
    cpmd[ln[nwise]]={}
    cpmd[ln[nwise]][samp]=ln[-1]
   else:
    cpmd[ln[nwise]][samp]=ln[-1]
 return cpmd
cpmd=CPMd(filelist,nwise)

def lgFC(cpmd,refd):
 tot=len(refd["0"])*len(refd["1"])
 fcd={}
 for com in cpmd:
  sfc=0
  for a in refd["0"]:   
   for b in refd["1"]:
    if a in cpmd[com] and b in cpmd[com]:
     key=b+"_"+a
     fc=float(cpmd[com][a])-float(cpmd[com][b]) 
     sfc+=fc
     if com not in fcd:
      fcd[com]={}
      fcd[com][key]=fc
     else:
      fcd[com][key]=fc
  if com in fcd:
   avfc=sfc/float(len(fcd[com])) 
   fcd[com]["avg"]=avfc
 return fcd
fcd=lgFC(cpmd,refd)

def npval(cpmd,refd):
 pvd={}
 for com in cpmd:
  pvd[com]=[[],[]]
  for a in refd["0"]:
   if a in cpmd[com]:
    pvd[com][0].append(float(cpmd[com][a]))
  for i in refd["1"]:
   if i in cpmd[com]:
    pvd[com][1].append(float(cpmd[com][i]))
 for com in pvd:
  a=np.array(pvd[com][0])
  b=np.array(pvd[com][1])
  (ts,pv)=stats.ttest_ind(a,b,equal_var=False)
  npv=-(math.log(pv+1))
  pvd[com].append(npv)
 return pvd
pvd=npval(cpmd,refd)

outfile=open("FCPV.csv","w")
row="key"
for i in range(nwise):
 row+=",guideRNA"+str(nwise-i)
row+=",log2FC,-log10pval"
outfile.write(row+"\r\n")
for key in pvd:
 row=key.strip("_")
 k=key.split("_")
 for i in k:
  row+=","+i
 if key in fcd and key in pvd:
  row+=","+str(fcd[key]["avg"])+","+str(pvd[key][-1])
  outfile.write(row+"\r\n")
 elif key in fcd and key not in pvd:
  row+=","+str(fcd[key]["avg"])+",nan" 
  outfile.write(row+"\r\n")
outfile.close()

