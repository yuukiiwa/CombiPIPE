import sys
fn=sys.argv[1]
ofn="sgRNAlevelGI.csv"
ofn2="genelevelGI.csv"

def addsgRNA(sgRNA,FC,sgDict):
 if sgRNA not in sgDict:
  sgDict[sgRNA]=[FC]
 else:
  sgDict[sgRNA].append(FC)
 return sgDict

def addCombo(ln,FC,comDict):
 combo=ln[0]
 cs=combo.split("+")
 revco=cs[1]+"+"+cs[0]
 if combo not in comDict and revco not in comDict:
  comDict[combo]=[FC]
 if combo in comDict:
  comDict[combo].append(FC)
 if revco in comDict:
  comDict[revco].append(FC)
 return comDict

def openFile(fn):
 sgDict,sgcomDict={},{}
 geDict,gecomDict={},{}
 file=open(fn,"r")
 header=file.readline().strip("\r\n").split(",")
 for ln in file:
  ln=ln.strip("\r\n").split(",")[1:]
  FC=float(ln[-2])
  #0->2; 1->3
  dummies=["1","2"]
  sgl,gl=ln,ln[1:]
  del sgl[1]
  if ln[2] in dummies and ln[1] not in dummies:
   addsgRNA(sgl[0].split("+")[0],FC,sgDict)
   addsgRNA(gl[0].split("+")[0],FC,geDict)
  if ln[1] in dummies and ln[2] not in dummies:
   addsgRNA(sgl[0].split("+")[1],FC,sgDict)
   addsgRNA(gl[0].split("+")[1],FC,geDict)
  if ln[1] not in dummies and ln[2] not in dummies:
   if ln[1] != ln[2]:
    addCombo(sgl,FC,sgcomDict)
    addCombo(gl,FC,gecomDict)
 for sgRNA in sgDict:
  sgDict[sgRNA].append(sum(sgDict[sgRNA])/len(sgDict[sgRNA]))
 for combo in sgcomDict:
  sgcomDict[combo].append(sum(sgcomDict[combo])/len(sgcomDict[combo]))
 for gene in geDict:
  geDict[gene].append(sum(geDict[gene])/len(geDict[gene]))
 for combo in gecomDict:
  gecomDict[combo].append(sum(gecomDict[combo])/len(gecomDict[combo]))
 return (sgDict,sgcomDict,geDict,gecomDict)
op=openFile(fn)
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
sgRNAgiDict=calcGI(sgDict,sgcomDict)
GENEgiDict=calcGI(geDict,gecomDict)

def outFile(giDict,ofn):
 outfile=open(ofn,"w")
 outfile.write("combo,expected,observed,GI"+"\r\n")
 for combo in giDict:
  row=combo
  for i in giDict[combo]:
   row+=","+str(i)
  outfile.write(row+"\r\n")
 outfile.close()
outFile(sgRNAgiDict,ofn)
outFile(GENEgiDict,ofn2)

def outDunnettFile(comDict,sgDict):
 for combo in comDict:
  t=combo.split("+")
  if t[0] != t[1]:
   outfile=open(combo+"_Dunnettin.csv","w")
   for fc in comDict[combo]:
    outfile.write(combo+",one_two,"+str(fc)+"\r\n")
   combo=combo.split("+")
   for sg in combo:
    for fc in sgDict[sg]:
     outfile.write(sg+",one,"+str(fc)+"\r\n")
   outfile.close()
outDunnettFile(gecomDict,geDict)
