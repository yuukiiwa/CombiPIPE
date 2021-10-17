import sys, os, math
from scipy import stats
import numpy as np
dir,sampleInfo,nwise=sys.argv[1],sys.argv[2],int(sys.argv[3])
def grabInfo(sampleInfo):
 file=open(sampleInfo,mode="r")
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
def CPMd(filelist):
 cpmd,named={},{}
 for fn in filelist:
  file=open("/"+dir+"/"+fn,"r")
  samp=fn.split("_")[1].split(".")[0]
  for y in range (8):
   file.readline()
  for ln in file:
   ln=ln.strip("\r\n").split(",")
   com=ln[2].strip("\ufeff")
   if com not in cpmd:
    cpmd[com]={}
    cpmd[com][samp]=ln[-1]
    named[com]=(ln[3],ln[4])
   else:
    cpmd[com][samp]=ln[-1]
 return (cpmd,named)
runCPMd=CPMd(filelist)
cpmd=runCPMd[0]
named=runCPMd[1]

def lgFC(cpmd,refd):
 tot=len(refd["0"])*len(refd["1"])
 fcd={}
 for com in cpmd:
  com=com.strip("\ufeff")
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
  avfc=sfc/float(len(fcd[com])) 
  fcd[com]["avg"]=avfc
 return fcd
fcd=lgFC(cpmd,refd)

def npval(cpmd,refd):
 pvd={}
 na = len(refd["0"])
 nb = len(refd["1"])
 for com in cpmd:
  pvd[com]=[[],[]]
  for a in refd["0"]:
   if a in cpmd[com]:
    pvd[com][0].append(float(cpmd[com][a]))
  for i in refd["1"]:
   if i in cpmd[com]:
    pvd[com][1].append(float(cpmd[com][i]))
 for com in pvd:
  if len(pvd[com][0]) == na and len(pvd[com][1]) == nb:
   a=np.array(pvd[com][0])
   b=np.array(pvd[com][1])
   (ts,pv)=stats.ttest_ind(a,b,equal_var=False)
   npv=-(math.log(pv))
   pvd[com].append(npv)
  else:
   pvd[com].append("") 
 return pvd
pvd=npval(cpmd,refd)

outfile=open("FCPV.csv","w")
row="key,sgRNAs,genes"
for i in range(nwise):
 row+=",guideRNA"+str(nwise-i)
row+=",log2FC,-log10pval"
outfile.write(row+"\r\n")
for key in pvd:
 row=key.strip("_")+","+named[key][0].strip("+")+","+named[key][1].strip("+")
 k=key.split("_")
 for i in k:
  row+=","+i
 row+=","+str(fcd[key]["avg"])+","+str(pvd[key][-1])
 outfile.write(row+"\r\n")
outfile.close()
