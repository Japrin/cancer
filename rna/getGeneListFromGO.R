#!/usr/bin/env Rscript

args <- commandArgs(T)
if(length(args)<2)
{
    cat("getGeneListFromGO.R <go.id> <out.file>\n")
    q()
}
go.id <- args[1]
out.file <- args[2]
#go.id  <- "GO:0005125"
#out.file <- "log"

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")
getEntrezFromGO <- function(term)
{
	if (require(org.Hs.eg.db)) 
	{
		xxGO <- AnnotationDbi::as.list(org.Hs.egGO2EG)
	}
	else
	{
		stop("Install org.Hs.eg.db package for retrieving gene lists from GO")
	}
   	eg <- unique(unlist(xxGO[term]))
    eg	
}

glist <- getEntrezFromGO(go.id)
o <- data.frame(GeneID=glist,GeneName=entrezToXXX(glist))
write.table(o,file = out.file,sep="\t",quote = F,col.names = F,row.names = F)
