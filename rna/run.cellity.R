#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--readCountFile", type="character", required=TRUE, help="read countt file ")
parser$add_argument("-d", "--designFile", type="character", required=TRUE, help="design matrix file;")
parser$add_argument("-a", "--statFile", type="character", required=TRUE, help="read stat file;")
parser$add_argument("-o", "--outDir", type="character", required=TRUE, help="output directory")
parser$add_argument("-e", "--ensembl", action="store_true", default=FALSE, help="geneID is ensembl id? [default %(default)s]")
parser$add_argument("-s", "--sampleID", type="character",default="SAMPLE", required=TRUE, help="sample id [default %(default)s]")
args <- parser$parse_args()
print(args)

library(cellity)
source("/Share/BP/zhenglt/02.pipeline/cancer/lib/scRNAToolKit.R")
source("/Share/BP/zhenglt/02.pipeline/cancer/rna/cellity.lib.R")
count.file <- args$readCountFile
design.file <- args$designFile
readStat.file <- args$statFile
out.dir <- args$outDir
sample.id <- args$sampleID
is.ensembl <- args$ensembl

#count.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/quantification/P0701.count.tab.gz"
#design.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/sample.design/P0701.above150K.sample.desc.txt"
#readStat.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/qc/cellity/OUT.stat/P0701.bamStat.forCellity.txt" 
#out.dir <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/qc/cellity/OUT.cellity"
#sample.id <- "P0701"

dir.create(out.dir,showWarnings = F,recursive = T)

myDesign<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
count.table <- read.table(count.file,header = T,row.names = 1,check.names = F,stringsAsFactors = F)
count.table <- count.table[,rownames(myDesign)]
if(!is.ensembl){
    cnames <- entrezToXXX(rownames(count.table),type = "ENSG")
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- rownames(count.table)[cnames.na]
    f.nodup <- !duplicated(cnames)
    count.nodup.table <- count.table[f.nodup,]
    rownames(count.nodup.table) <- cnames[f.nodup]
}else{
    count.nodup.table <- count.table
    rownames(count.nodup.table) <- gsub("\\.\\d+$","",rownames(count.nodup.table),perl=T)
}


readStat.table <- read.csv(readStat.file,header=T,check.names=F,stringsAsFactors = F)
rownames(readStat.table) <- readStat.table[,1]
readStat.table <- readStat.table[rownames(myDesign),]

#### 
count.nodup.table.nm <- normalise_by_factor(count.nodup.table, colSums(count.nodup.table))
print(head(count.nodup.table.nm[,1:4]))
print(head(readStat.table[,1:4]))
## some bug with this version
##sample_features <- extract_features(count.nodup.table.nm, readStat.table,organism="human",output_dir = out.dir,prefix = sample.id)
sample_features <- my_extract_features(count.nodup.table.nm, readStat.table,organism="human",output_dir = out.dir,prefix = sample.id)

sample_features_all <- sample_features[[1]]
pdf(sprintf("%s/%s.cellity.pdf",out.dir,sample.id),width = 12,height = 8)
sample_qual_pca <- assess_cell_quality_PCA(sample_features_all,file=sprintf("%s/%s.cellity.pdf",out.dir,sample.id))
dev.off()
write.table(sample_qual_pca,file=sprintf("%s/%s.cellity.txt",out.dir,sample.id),quote = F,sep = "\t",row.names = F)

