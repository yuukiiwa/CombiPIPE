import sys
from itertools import permutations

fn,nwise,dummysg=sys.argv[1],int(sys.argv[2]),sys.argv[3]
dummysg=dummysg.split(",")
ofn="sgRNAlevelGI.csv"
ofn2="genelevelGI.csv"

def addsgRNA(sgRNA,FC,sgDict):
 if sgRNA not in sgDict:
  sgDict[sgRNA]=[FC]
 else:
  sgDict[sgRNA].append(FC)
 return sgDict

def addCombo(com,FC,comDict):
 combos=list(permutations(com))
 combos=["+".join(combo) for combo in combos]
 if not any(combo in comDict for combo in combos): 
  comDict[combos[0]]=[FC]
 else:
  for combo in combos:
   if combo in comDict:
    comDict[combo].append(FC)
 return comDict

def openFile(fn,niwse,dummies):
 sgDict,geDict,sgcomDict,gecomDict={},{},{},{}
 file=open(fn,"r")
 header=file.readline().strip("\r\n").split(",")
 for ln in file:
  ln=ln.strip("\r\n").split(",")[1:]
  FC=float(ln[-2])
  #0->2; 1->3
  sgcombo,gecombo,dcont=ln[0].split("+"),ln[1].split("+"),0
  for d in dummies:
   dcont+=sgcombo.count(d)
  if dcont == nwise-1:
   for sg in sgcombo:
    if sg not in dummies:
     sgRNA=sg
   addsgRNA(sgRNA,FC,sgDict)
   for ge in gecombo:
    if ge not in dummies:
     gene=ge
   addsgRNA(gene,FC,geDict)
  else:
   if len([x for x in sgcombo if sgcombo.count(x) > 1]) == 0:
    addCombo(sgcombo,FC,sgcomDict)
   if len([x for x in gecombo if gecombo.count(x) > 1]) == 0:
    addCombo(gecombo,FC,gecomDict)
 for sgRNA in sgDict:
  sgDict[sgRNA].append(sum(sgDict[sgRNA])/len(sgDict[sgRNA]))
 for sgcombo in sgcomDict:
  sgcomDict[sgcombo].append(sum(sgcomDict[sgcombo])/len(sgcomDict[sgcombo]))
 for gene in geDict:
  geDict[gene].append(sum(geDict[gene])/len(geDict[gene]))
 for gecombo in gecomDict:
  gecomDict[gecombo].append(sum(gecomDict[gecombo])/len(gecomDict[gecombo]))
 return (sgDict,sgcomDict,geDict,gecomDict)
op=openFile(fn,nwise,dummysg)
sgDict,sgcomDict,geDict,gecomDict=op[0],op[1],op[2],op[3]  #the last item is the mean FC

def calcGI(sgDict,comDict):
 giDict={}
 for combo in comDict:
  sgs=combo.split("+")
  obs=comDict[combo][-1]
  exp=0
  for sg in sgs:
   if sg in sgDict:
    exp+=sgDict[sg][-1]
  GI=obs-exp
  giDict[combo]=(exp,obs,GI)
 return giDict
sgGIDict=calcGI(sgDict,sgcomDict)
geGIDict=calcGI(geDict,gecomDict)

def outFile(giDict,ofn):
 outfile=open(ofn,"w")
 outfile.write("gene combo,expected,observed,GI"+"\r\n")
 for combo in giDict:
  row=combo
  for i in giDict[combo]:
   row+=","+str(i)
  outfile.write(row+"\r\n")
 outfile.close()
outFile(sgGIDict,ofn)
outFile(geGIDict,ofn2)

def outDunnettFile(comDict,sgDict):
 for combo in comDict:
  outfile=open(combo+"_Dunnettin.csv","w")
  for fc in comDict[combo]:
   outfile.write(combo+",one_two,"+str(fc)+"\r\n")
  combo=combo.split("_")
  for sg in combo:
   if sg in sgDict:
    for fc in sgDict[sg]:
     outfile.write(sg+",one,"+str(fc)+"\r\n")
  outfile.close()
outDunnettFile(gecomDict,geDict)
