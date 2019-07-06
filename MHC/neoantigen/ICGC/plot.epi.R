#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<2)
{
    cat("plot.epi.R <in.file> <out.prefix>\n")
    q()
}


library(beeswarm)
library(reshape2)
#library(ggplot2)

in.file <- args[1]
out.prefix <- args[2]
#in.file <- "HLA-A0201.9.netMHC.out.filtered.pop.summary.tab"
#out.prefix <- "HLA-A0201.9.netMHC.out.filtered.pop.summary.tab"
print(in.file)
in.table <- read.table(in.file,header=T,check.names=F,stringsAsFactors = F)
in.table.melt <- melt(in.table,id.vars=c("CancerType","SampleID"),variable.name="CountType",value.name="Count")
in.table$PercentOfEpi  <- in.table$Epi*100 / in.table$MutBurden

## sort by median
g <- aggregate(in.table$Epi, by=list(CancerType=in.table$CancerType), FUN=median)
g.percent <- aggregate(in.table$PercentOfEpi, by=list(CancerType=in.table$CancerType), FUN=median)
g.nsample <- aggregate(in.table$PercentOfEpi, by=list(CancerType=in.table$CancerType), FUN=length)
g.a <- merge(g.percent,g.nsample,by="CancerType")
g <- merge(g,g.a,by="CancerType")
names(g) <- c("CancerType","EpiMedian","EpiPercent","NSample")
g <- g[order(g$EpiMedian,decreasing = T),]
print(g)
in.table.melt$CancerType <- factor(x=in.table.melt$CancerType ,levels=g$CancerType,ordered=T)
in.table$CancerType <- factor(x=in.table$CancerType ,levels=g$CancerType,ordered=T)
head(in.table.melt)

## mydotplot
myDotPlot <- function(dat,midpoints)
{
  dat <- dat[order(dat$Count,decreasing = F),]
  print(head(dat))

  g1 <- aggregate(dat$CancerType,by=list(CancerType=dat$CancerType),length)
  g2 <- aggregate(dat$Count,by=list(CancerType=dat$CancerType),median)
  g <- merge(g1,g2,by="CancerType")
  names(g) <- c("CancerType","NumSample","NumEpi")
  for(i in 1:nrow(g))
  {
    j <- as.numeric(g[i,1])
    n <- g[i,2]
    m <- g[i,3]
    segments(midpoints[j]-0.5,m,midpoints[j]+0.5,m)
    dat.p <- subset(dat,CancerType==g[i,1])
    y.p <- dat.p$Count
    x.p <- midpoints[j]-0.5+(0:(n-1))/(n-1)
    points(x.p,y.p,pch=20,cex=0.02)
    cat(paste0("CancerType: ",g[i,1],",n=",n,"\n"))
    print(y.p)
    
    
  }
}

#png(paste0(out.prefix,".png"),width=2000,height=800)
pdf(paste0(out.prefix,".pdf"),width=20,height=8)
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE),heights = c(1.5,3,3,1))
par(mar=c(3,10,1,3),cex.axis=2.5,cex.main=2,cex.lab=2.8)
## barplot
dat.plot <- t(g$EpiPercent)
colnames(dat.plot) <- g$CancerType
xx<-barplot(dat.plot,ylab="neoantigen\nfraction(%)\n",xaxt="n",main="",cex.lab=2)
#text(xx,par("usr")[3]-1, srt = 45, adj = 1, labels = paste0(g$CancerType,"\n(",g$NSample,")"), xpd = TRUE,cex=2.8)
## mydotplot 1
par(mar=c(2.5,10,1,3),cex.axis=2.5,cex.main=2,cex.lab=2.8)
dat.plot <- rep(50,length(dat.plot))
names(dat.plot) <- g$CancerType
xx<-barplot(dat.plot,ylab="#neoantigen\n",border=NA,col="white", ,xaxt="n")
#text(xx,par("usr")[3]-10, srt = 45, adj = 1, labels = paste0(g$CancerType,"\n(",g$NSample,")"), xpd = TRUE,cex=2.8)
myDotPlot(subset(in.table.melt,CountType=="Epi"),xx)
## mydotplot 2
par(mar=c(2.5,10,1,3),cex.axis=2.5,cex.main=2,cex.lab=2.8)
dat.plot <- rep(500,length(dat.plot))
names(dat.plot) <- g$CancerType
xx<-barplot(dat.plot,ylab="#mutation\n",border=NA,col="white", ,xaxt="n")
#text(xx,par("usr")[3]-100, srt = 45, adj = 1, labels = paste0(g$CancerType,"\n(",g$NSample,")"), xpd = TRUE,cex=2.8)
myDotPlot(subset(in.table.melt,CountType=="MutBurden"),xx)
xx<-barplot(dat.plot,plot = FALSE,ylab="",border=NA,col="white", ,xaxt="n")
text(xx,par("usr")[3]-50, srt = 45, adj = 1, labels = paste0(g$CancerType,"\n(",g$NSample,")"), xpd = NA,cex=2.8)
dev.off()











