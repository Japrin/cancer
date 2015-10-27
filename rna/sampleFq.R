#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ShortRead"))

parser <- ArgumentParser()
parser$add_argument("-n", "--nreads", type="double", default=1e6, help="number of reads to sample [default %(default)s]")
parser$add_argument("-s", "--seed", type="integer", default=11, help="seed of random number [default %(default)s]")
parser$add_argument("infile", help="input file")
parser$add_argument("outfile", help="output file")
args <- parser$parse_args()

infile <- args$infile
outfile <- args$outfile
n <- args$n
s <- args$seed

set.seed(s)
sampler <- FastqSampler(infile,n=n)
fq <- yield(sampler)
r=writeFastq(fq,outfile,compress=T)
print(r)
close(sampler)

