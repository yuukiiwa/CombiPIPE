import sys, os
dir,fcgi=sys.argv[1],open(sys.argv[2],"r")
outfile=open("genlevelGIdunnettPV.csv","w")
filelist=os.listdir(dir)
pdict={}
for fn in filelist:
 file=open("/"+dir+"/"+fn,"r")
 file.readline()
 fn=fn[7:-14]
 pval=file.readline().split(",")[-1]
 pdict[fn]=pval.strip("\n")
print(pdict)
outfile.write(fcgi.readline().strip("\r\n")+",Dunnett p-val"+"\r\n")
for ln in fcgi:
 c=ln.split(",")[0]
 if c in pdict:
  outfile.write(ln.strip("\r\n")+","+pdict[c]+"\r\n")
outfile.close()
