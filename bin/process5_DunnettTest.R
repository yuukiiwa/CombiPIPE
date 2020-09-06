args = commandArgs(trailingOnly=TRUE)
path<-args[1]
files <- grep(list.files(path), pattern='genelevelGI.csv', inv=T, value=T)

library(DescTools)
library(ggplot2)
for(i in 1:length(files)){
  if (files[i] != "sgRNAlevelGI.csv"){
    fullpath<-paste(path,files[i],sep='/') 
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
    write.csv(dunnett, file = reportname)
    ####ploting
    print(files[i])
    log2FC<-df$fold_change
    Combinations<-df$combination
    Type<-df$combinations
    p <- ggplot(df, aes(x=combination, y=fold_change, color=combinations)) +
         geom_violin() +
         geom_point(size=4, alpha=0.7, position=position_jitter(w=0.1, h=0)) +
         stat_summary(fun.y=mean, geom="point", shape=23, color="black", size=4) +         
         stat_summary(fun.ymin=function(x)(mean(x)-sd(x)), 
                      fun.ymax=function(x)(mean(x)+sd(x)),
                      geom="errorbar", width=0.1) +
                      theme_bw()
    imagename <- paste("DunnettPlot",files[i],".png",sep="") 
    ggsave(imagename,plot=p)
  }
}
