#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--in",required=T,help="input file list,first line is header with filename, sample")
parser$add_argument("-o", "--out",required=T,help="output prefix")
parser$add_argument("-t", "--txfile",default="/DBS/DB_temp/zhangLab/ensemble/mybuild/kallisto/b20170519/rel88/ID.mapping.txt",
                    help="tx to gene file, no header with first 3 column txID, geneID and geneSymbol [default %(default)s]")
parser$add_argument("-a", "--txOut", action="store_true", default=FALSE, help="output transcript level expression [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose output [default %(default)s]")
args <- parser$parse_args()
args.txOut <- args$txOut

print(args)

in.file <- args[["in"]]
out.prefix <- args[["out"]]
tx.file <- args[["txfile"]]

dir.create(dirname(out.prefix),showWarnings = F,recursive = T)
#in.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/quantification/test.kallisto.list"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/quantification/test.tximport.out"
#tx.file <- "/DBS/DB_temp/zhangLab/ensemble/mybuild/kallisto/ID.mapping.txt"

library(tximport)
library(readr)
#options(stringsAsFactors = F)

sample.desc <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F)
files <- sample.desc$filename
names(files) <- sample.desc$sample
tx.desc <- read.table(tx.file,header = F,check.names = F,stringsAsFactors = F,sep="\t")[,1:3]
rownames(tx.desc) <- tx.desc[,1]
colnames(tx.desc) <- c("TXNAME","GENEID","GENESYMBOL")

gene.desc <- unique(tx.desc[,c(2,3)])
rownames(gene.desc) <- gene.desc[,1]

#### tximport_1.4.0
## some bug (https://support.bioconductor.org/p/87553/)
#txi <- tximport(files, type = "kallisto", tx2gene = tx.desc, importer=read_tsv)
#txi <- tximport(files, type = "kallisto", tx2gene = tx.desc, importer=read.delim)
nn <- nrow(read.delim(files[1]))
my.read.tsv <- function(x) readr::read_tsv(x, guess_max=nn)
txi <- tximport(files, type = "kallisto", tx2gene = tx.desc, importer=my.read.tsv,txOut=args.txOut)
####
##txi <- tximport(files, type = "kallisto", tx2gene = tx.desc, reader = read_tsv)
names(txi)

for(cp in c("abundance","counts","length"))
{
    ##out.df <- data.frame(geneID=rownames(txi[[cp]]))
    out.df <- data.frame(geneID=rownames(txi[[cp]]),stringsAsFactors = F)
    out.df$geneSymbol <- gene.desc[ out.df$geneID,"GENESYMBOL"]
    out.df <- cbind(out.df,txi[[cp]])
    conn <- gzfile(sprintf("%s.kallisto.tximport.%s.txOut%s.txt.gz",out.prefix,cp,if(args.txOut) "T" else "F"),open = "w")
    write.table(out.df,conn,quote = F,sep = "\t",row.names = F)
    close(conn)
}

save(txi,file = sprintf("%s.kallisto.tximport.RData",out.prefix))

cat("run.txImport.R done\n")

