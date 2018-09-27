#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("R.utils"))

parser <- ArgumentParser()
parser$add_argument("-i", "--infile",required=T,help="input file")
parser$add_argument("-o", "--outprefix",required=T,help="output prefix")
parser$add_argument("-d", "--samplefile",required=T,help="sample description file")
parser$add_argument("-p", "--prop",default=0.1,help="mt exp threshold [default %(default)s]")
parser$add_argument("-s", "--sampleID",default="TSAMPLE",help="sample id [default %(default)s]")
parser$add_argument("-t", "--txfile",default="/DBS/DB_temp/zhangLab/ensemble/mybuild/kallisto/ID.mapping.MT.txt",
                    help="tx to gene file, no header with first 3 column txID, geneID and geneSymbol [default %(default)s]")
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose output [default %(default)s]")
args <- parser$parse_args()

print(args)

in.file <- args[["infile"]]
out.prefix <- args[["outprefix"]]
tx.file <- args[["txfile"]]
d.file <- args[["samplefile"]]
propT <- args[["prop"]]
sample.id <- args[["sampleID"]]

#in.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/quantification/test.tximport.out.kallisto.tximport.abundance.txt.gz"
#out.prefix <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/quantification/test.mt"
#tx.file <- "/DBS/DB_temp/zhangLab/ensemble/mybuild/kallisto/ID.mapping.MT.txt"
#d.file <- "/WPSnew/zhenglt/work/proj_zl/integrate0924/quantification/test.kallisto.list"
#propT <- 0.1
#sample.id <- "TEST"

#options(stringsAsFactors = F)
### exp file
if(grepl(".txt.gz",in.file,perl = T)){
    in.table <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F)
    rownames(in.table) <- in.table[,1]
    g.GNAME <- in.table[,2]
    names(g.GNAME) <- rownames(in.table)
    in.table <- in.table[,c(-1,-2)]
}else{
    lenv <- loadToEnv(in.file)
    in.table <- lenv$Y
    g.GNAME <- lenv$g.GNAME
}

### sample file
myDesign <- read.table(d.file,header = T,check.names = F,stringsAsFactors = F)
rownames(myDesign) <- myDesign$sample
f.sample <- intersect(colnames(in.table),myDesign$sample)
in.table <- in.table[,f.sample]
myDesign <- myDesign[f.sample,]
### MT gene file
tx.desc <- read.table(tx.file,header = F,check.names = F,stringsAsFactors = F,sep="\t")[,1:3]
rownames(tx.desc) <- tx.desc[,1]
colnames(tx.desc) <- c("TXNAME","GENEID","GENESYMBOL")
### calcualte the proportion
##mt.prop <- apply(in.table[tx.desc$GENEID,],2,function(x){ sum(x)/1e6 })
mt.prop <- apply(in.table[tx.desc$GENEID,],2,function(x){ sum(x) })
e.total <- apply(in.table,2,sum)
mt.prop <- mt.prop/e.total

myDesign$MT.prop <- mt.prop
write.table(myDesign,file = sprintf("%s.Extended.txt",out.prefix),row.names = F,quote = F,sep = "\t")
write.table(subset(myDesign,MT.prop<=propT),file = sprintf("%s.flt.txt",out.prefix),row.names = F,quote = F,sep = "\t")
pdf(sprintf("%s.propDist.pdf",out.prefix),width = 6,height = 6)
hist(myDesign$MT.prop,breaks = 50,xlab="Prop by Mitochondria Genes",main=sample.id)
abline(v=propT,lty=2,lw=1.5)
dev.off()

cat("qc.filterByMTProp.R done.\n")

