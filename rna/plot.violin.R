#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-i", "--in",required=T)
parser$add_argument("-o", "--out",required=T)
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE, help="verbose output [default %(default)s]")
parser$add_argument("-s", "--sample", default="SAMPLE", help="sample id [default %(default)s]")
parser$add_argument("-f", "--genelist",type="character", help="gene list file, format: geneID[tab]geneSymbol [default %(default)s]")
args <- parser$parse_args()

in.bam.file <- args$in
out.prefix <- args$out

