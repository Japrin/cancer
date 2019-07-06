#!/usr/bin/env Rscript

args <- commandArgs(T)

if(length(args)<2)
{
	cat("usage: binder.count.R <in.file> <out.prefix>\n")
	q()
}

in.file <- args[1]
out.prefix <- args[2]

#library(beeswarm)
library(reshape2)
#library(ggplot2)

#in.file <- "/WPS1/zhenglt/work/MHC/TCGA/OUT.netMHC-3.4/HLA-A02:01.9/HLA-A02:01.9.netMHC.out.filtered.stat.txt"
#out.prefix <- "aaaa.binder"

in.table <- read.table(in.file,header=T,check.names=F,stringsAsFactors = F)
in.table.melt <- melt(in.table,id.vars=c("CancerType","Sample"),variable.name="CountType",value.name="Count")
in.table$PercentOfBinder  <- in.table$NumOfBinder / in.table$NumOfMissense

g <- aggregate(in.table$NumOfMissense, by=list(CancerType=in.table$CancerType), FUN=median)
g <- g[order(g$x,decreasing = T),]
in.table.melt$CancerType <- factor(x=in.table.melt$CancerType ,levels=g$CancerType,ordered=T)
in.table$CancerType <- factor(x=in.table$CancerType ,levels=g$CancerType,ordered=T)

png(paste0(out.prefix,".png"),width=1500,height=800)
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE))
par(mar=c(5,4,4,2),cex.axis=2,cex.main=2)
boxplot(Count~CancerType,data=in.table.melt[in.table.melt$CountType=="NumOfMissense",], outline = TRUE,ylim=c(0,1200),main="Number Of Missense")
#par(mar=c(5,4,4,2))
boxplot(Count~CancerType,data=in.table.melt[in.table.melt$CountType=="NumOfBinder",], outline = TRUE,ylim=c(0,150),main="Number Of Binder")
#par(mar=c(5,4,4,2))
boxplot(PercentOfBinder~CancerType,data=in.table,main="Binder Percentage",ylim=c(0,0.5))
dev.off()

