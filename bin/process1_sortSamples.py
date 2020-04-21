import regex, sys
fastq,sampinfo,identifier_pattern,dict=sys.argv[1],sys.argv[2],sys.argv[3],{}
sfile=open(sampinfo,"r")
for ln in sfile:
 ln=ln.strip("\r\n").split(",")
 dict[ln[0]]=ln[1]

def oneMismatch(id,seq):
 restring='('+id+'){e<=1}'
 m=regex.findall(restring,seq,overlapped=False)
 return m

file=open(fastq,"r")
line,lines=[],[]
for ln in file:
 if ln.startswith(identifier_pattern):
  if len(line) != 0:
   lines.append(line)
   line=[]
  line.append(ln)
 else:
  line.append(ln)
lines.append(line)

for id in dict:
 n=len(id)+2
 ofn=id+"_"+dict[id]+".fastq"
 outfile=open(ofn,"w")
 for line in lines:
  m=oneMismatch(id,line[1][:n])
  if len(m) > 0:
   for ln in line:
    outfile.write(ln)
 outfile.close()
