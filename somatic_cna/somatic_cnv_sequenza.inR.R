#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--inFile", type="character", required=TRUE, help="input file")
parser$add_argument("-o", "--outPrefix", type="character", required=TRUE, help="outPrefix")
parser$add_argument("-s", "--sampleID", type="character",default="TEST", help="sampleID [default %(default)s]")
parser$add_argument("-x", "--sex", type="character",default="F", help="sex, M or F [default %(default)s]")

args <- parser$parse_args()
print(args)

############## tune parametrs  ########

in.file <- args$inFile
out.prefix <- args$outPrefix
sampleID <- args$sampleID
opt.sex <- args$sex

dir.create(dirname(out.prefix),F,T)

library("sequenza")
library("tictoc")

tic("sequenza.extract")
se.list <- sequenza.extract(in.file,
                         parallel=8,
                         assembly="hg38",
                         verbose = T)
toc()

tic("sequenza.fit")
CP.list <- sequenza.fit(se.list,
                   female=(opt.sex=="F"),
                   XY = c(X = "chrX", Y = "chrY"),
                   mc.cores=8)
toc()

tic("sequenza.results")
sequenza.results(sequenza.extract = se.list,
                 cp.table = CP.list,
                 female=(opt.sex=="F"),
                 XY = c(X = "chrX", Y = "chrY"),
                 sample.id = sampleID, out.dir=dirname(out.prefix))
toc()


