#!/usr/local/bin/Rscript

args<-commandArgs(T)
if(length(args)<1)
{
	cat("usage: CNVseq.print.CNV.R <infile>\n")
	q()
}

library(cnv)

infile<-args[1]

data <- read.delim(infile)
cnv.print(data)
