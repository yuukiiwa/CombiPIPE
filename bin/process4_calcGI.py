import sys
fn=sys.argv[1]
ofn="genelevelGI.csv"

def addsgRNA(sgRNA,FC,sgDict):
 if sgRNA not in sgDict:
  sgDict[sgRNA]=[FC]
 else:
  sgDict[sgRNA].append(FC)
 return sgDict

def addCombo(ln,FC,comDict):
 combo=ln[0]+"-"+ln[1]
 revco=ln[1]+"-"+ln[0]
 if combo not in comDict and revco not in comDict:
  comDict[combo]=[FC]
 if combo in comDict:
  comDict[combo].append(FC)
 if revco in comDict:
  comDict[revco].append(FC)
 return comDict

def openFile(fn):
 sgDict, comDict={},{}
 file=open(fn,"r")
 header=file.readline().strip("\r\n").split(",")
 for ln in file:
  ln=ln.strip("\r\n").split(",")[1:]
  FC=float(ln[-2])
  #0->2; 1->3
  dummies=["1","2"]
  if ln[0] in dummies and ln[1] not in dummies:
   addsgRNA(ln[1],FC,sgDict)
  if ln[1] in dummies and ln[0] not in dummies:
   addsgRNA(ln[0],FC,sgDict)
  if ln[1] not in dummies and ln[0] not in dummies:
   if ln[1] != ln[0]:
    addCombo(ln,FC,comDict)
 for sgRNA in sgDict:
  sgDict[sgRNA].append(sum(sgDict[sgRNA])/len(sgDict[sgRNA]))
 for combo in comDict:
  comDict[combo].append(sum(comDict[combo])/len(comDict[combo]))
 return (sgDict,comDict)
op=openFile(fn)
sgDict,comDict=op[0],op[1]  #the last item is the mean FC
print(comDict)

def calcGI(sgDict,comDict):
 giDict={}
 for combo in comDict:
  sgs=combo.split("-")
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
  combo=combo.split("-")
  for sg in combo:
   for fc in sgDict[sg]:
    outfile.write(sg+",one,"+str(fc)+"\r\n")
  outfile.close()
outDunnettFile(comDict,sgDict)
