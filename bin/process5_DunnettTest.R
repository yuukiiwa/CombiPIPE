args = commandArgs(trailingOnly=TRUE)
path<-args[1]
files <- grep(list.files(path), pattern='FCPVGI.csv', inv=T, value=T)

library(DescTools)
for(i in 1:length(files)){
  fullpath<-paste(path,files[1],sep='/') 
  counts <- read.table(fullpath,sep=",")
  colnames(counts)<- c("combination","combinations","fold_change")
  df <- counts
  df$combinations <- as.factor(df$combinations)
  mod1= lm(fold_change ~ combinations, data= df)
  summary(mod1)
  anova(mod1)
  confint(mod1)
  btwn <- aov(fold_change ~ combinations, data=df)
  summary(btwn)
  res<- DunnettTest(fold_change ~ combinations, data=df, control="one_two")
  dunnett= as.data.frame(res[["one_two"]])
  reportname= paste("Dunnett",files[i],sep="")
  print(reportname)
  write.csv(dunnett, file = reportname)
}

