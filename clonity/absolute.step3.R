#!/usr/bin/env Rscript

library(ABSOLUTE)
obj.name <- "ESCC.S.sample17"
calls.path = "/Share/BP/zhenglt/shidao/ExomeSeq/analysis/tumor-tissue/absolute.all/ESCC.S.sample17.PP-calls_tab.reviewed.txt"
modes.path = "/Share/BP/zhenglt/shidao/ExomeSeq/analysis/tumor-tissue/absolute.all/ESCC.S.sample17.PP-modes.data.RData"
output.path = "/Share/BP/zhenglt/shidao/ExomeSeq/analysis/tumor-tissue/absolute.all"

ExtractReviewedResults(calls.path, "japrin", modes.path, output.path, obj.name, "total", verbose=TRUE)
