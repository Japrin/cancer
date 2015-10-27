#!/usr/bin/env Rscript

library(ABSOLUTE)
obj.name <- "ESCC.S.sample17"
rda.files <- read.table("absolute.rdata.files",header=F,stringsAsFactors=F,check.names=F)
names(rda.files) <- c("sampleID","rda")
outDir <- "/Share/BP/zhenglt/shidao/ExomeSeq/analysis/tumor-tissue/absolute.all"

CreateReviewObject(obj.name, rda.files$rda, outDir, "total", verbose=TRUE)
