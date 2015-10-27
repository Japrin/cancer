#!/usr/bin/env Rscript

cyclebase.rdata <- "./cycleBase.human.RData"
lname <- load(cyclebase.rdata)
require("hom.Mm.inp.db")
aa<-inpIDMapper(ids=(dataCB$gene_systematic_name),"HOMSA", "MUSMU",srcIDType="ENSEMBL",destIDType="ENSEMBL")
dataCB$mouse.ENSG <- aa[as.character(dataCB$gene_systematic_name) ]

aa<-inpIDMapper(ids=(dataCB$gene_systematic_name),"HOMSA", "MUSMU",srcIDType="ENSEMBL",destIDType="SYMBOL")
dataCB$mouse.geneSymbol <- aa[as.character(dataCB$gene_systematic_name) ]

bb<-inpIDMapper(ids=(dataCB$gene_systematic_name),"HOMSA", "MUSMU",srcIDType="ENSEMBL",destIDType="EG")
dataCB$mouse.entrez <- bb[as.character(dataCB$gene_systematic_name) ]

save(dataCB,file = cyclebase.rdata)
