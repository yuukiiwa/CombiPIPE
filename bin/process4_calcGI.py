import sys
from itertools import permutations

fn,nwise,dummysg=sys.argv[1],int(sys.argv[2]),sys.argv[3]
dummysg=dummysg.split(",")
ofn="genelevelGI.csv"

def addsgRNA(sgRNA,FC,sgDict):
 if sgRNA not in sgDict:
  sgDict[sgRNA]=[FC]
 else:
  sgDict[sgRNA].append(FC)
 return sgDict

def addCombo(com,FC,comDict):
 combos=list(permutations(com))
 combos=["_".join(combo) for combo in combos]
 if not any(combo in comDict for combo in combos): 
  comDict[combos[0]]=[FC]
 else:
  for combo in combos:
   if combo in comDict:
    comDict[combo].append(FC)
 return comDict

def openFile(fn,niwse,dummies):
 sgDict, comDict={},{}
 file=open(fn,"r")
 header=file.readline().strip("\r\n").split(",")
 for ln in file:
  ln=ln.strip("\r\n").split(",")[1:]
  FC=float(ln[-2])
  #0->2; 1->3
  combo,dcont=ln[:nwise],0
  for d in dummies:
   dcont+=combo.count(d)
  if dcont == nwise-1:
   for sg in combo:
    if sg not in dummies:
     sgRNA=sg
   addsgRNA(sgRNA,FC,sgDict)
  else:
   if len([x for x in combo if combo.count(x) > 1]) == 0:
    addCombo(combo,FC,comDict)
 for sgRNA in sgDict:
  sgDict[sgRNA].append(sum(sgDict[sgRNA])/len(sgDict[sgRNA]))
 for combo in comDict:
  comDict[combo].append(sum(comDict[combo])/len(comDict[combo]))
 return (sgDict,comDict)
op=openFile(fn,nwise,dummysg)
sgDict,comDict=op[0],op[1]  #the last item is the mean FC

def calcGI(sgDict,comDict):
 giDict={}
 for combo in comDict:
  sgs=combo.split("_")
  obs=comDict[combo][-1]
  exp=0
  for sg in sgs:
   if sg in sgDict:
    exp+=sgDict[sg][-1]
  GI=obs-exp
  giDict[combo]=(exp,obs,GI)
 return giDict
giDict=calcGI(sgDict,comDict)

def outFile(giDict,ofn):
 outfile=open(ofn,"w")
 outfile.write("gene combo,expected,observed,GI"+"\r\n")
 for combo in giDict:
  row=combo
  for i in giDict[combo]:
   row+=","+str(i)
  outfile.write(row+"\r\n")
 outfile.close()
outFile(giDict,ofn)

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
outDunnettFile(comDict,sgDict)
