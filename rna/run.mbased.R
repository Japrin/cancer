#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("-t", "--type", type="integer", default=2, help="2 samples OR 1 sample analysis [default %(default)s]",metavar="number")
#parser$add_argument("-n", "--add_numbers", action="store_true", default=FALSE, help="Print line number at the beginning of each line [default]")
parser$add_argument("infile", nargs=1, help="Infile")
args <- parser$parse_args()
infile <- args$infile

if( file.access(infile) == -1) {
	stop(sprintf("Specified file ( %s ) does not exist", infile))
}

library(MBASED)

###### function definition ######
summarizeASEResults_2s <- function(MBASEDOutput) {
    geneOutputDF <- data.frame(
    	majorAlleleFrequencyDifference=assays(MBASEDOutput)$majorAlleleFrequencyDifference[,1],
    	pValueASE=assays(MBASEDOutput)$pValueASE[,1],
    	pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
    )
    geneOutputDF$pValueASE.adj <- p.adjust(geneOutputDF$pValueASE,method="BH")
    geneOutputDF$gene<-row.names(geneOutputDF)

    lociOutputGR <- rowData(exptData(MBASEDOutput)$locusSpecificResults)
    lociOutputGR$allele1IsMajor <- assays(exptData(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
    lociOutputGR$MAFDifference <- assays(exptData(MBASEDOutput)$locusSpecificResults)$MAFDifference[,1]
    lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID))) 
    return(
    	list(
    		geneOutput=geneOutputDF,
    		locusOutput=lociOutputList
    	)
    )
}
summarizeASEResults_1s <- function(MBASEDOutput) {
	geneOutputDF <- data.frame(
		majorAlleleFrequency=assays(MBASEDOutput)$majorAlleleFrequency[,1],
        	pValueASE=assays(MBASEDOutput)$pValueASE[,1],  
	   	pValueHeterogeneity=assays(MBASEDOutput)$pValueHeterogeneity[,1]
	)   
    	geneOutputDF$pValueASE.adj <- p.adjust(geneOutputDF$pValueASE,method="BH")
    	geneOutputDF$gene<-row.names(geneOutputDF)

	lociOutputGR <- rowData(exptData(MBASEDOutput)$locusSpecificResults)
	lociOutputGR$allele1IsMajor <- assays(exptData(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
	lociOutputGR$MAF <- assays(exptData(MBASEDOutput)$locusSpecificResults)$MAF[,1]
	lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID))) 
	return(
		list(
			geneOutput=geneOutputDF,
			locusOutput=lociOutputList
		)
	)
}


inTable<-read.table(infile,header=T,check.names=F,stringsAsFactors=F)
#colnames(inTable)<-c("seqnames","ranges","allele1","allele2","aseID","lociAllele1Counts_normal","lociAllele2Counts_normal","lociAllele1Counts_tumor","lociAllele2Counts_tumor","row_id")
inTable_rowData<-GRanges(seqnames=inTable$seqnames,ranges=IRanges(start=inTable$ranges, width=1),aseID=inTable$aseID,allele1=inTable$allele1,allele2=inTable$allele2)
names(inTable_rowData)<-inTable$row_id

## create input SummarizedExperiment object
if(args$type==1)
{
    mySEData <- SummarizedExperiment(
	    assays=list(
		    lociAllele1Counts=matrix(
			    inTable[,"lociAllele1Counts"],
			    ncol=1,
			    dimnames=list(
				    names(inTable_rowData), 
				    'tumor'
			    )
		    ),
		    lociAllele2Counts=matrix(
			    inTable[,"lociAllele2Counts"],
			    ncol=1, 
			    dimnames=list(
				    names(inTable_rowData), 
				    'tumor'
			    )
		    )
	    ),
	    rowData=inTable_rowData
    )
    ASEresults <- runMBASED(
    	ASESummarizedExperiment=mySEData,
    	isPhased=FALSE,
    	numSim=10^6,
    	BPPARAM = MulticoreParam()
    )
    ASEresult_1s_sum <- summarizeASEResults_1s(ASEresults)
    f  <- with(ASEresult_1s_sum$geneOutput, abs(majorAlleleFrequency-0.5) > 0.2 & pValueASE.adj < 0.05)
    f2 <- with(ASEresult_1s_sum$geneOutput, abs(majorAlleleFrequency-0.5) > 0.2 & pValueASE.adj < 0.05 & !is.na(pValueHeterogeneity) & pValueHeterogeneity <0.01 )
    write.table(ASEresult_1s_sum$geneOutput,paste(infile,".ASE.1s.Gene",sep=""),sep="\t",quote=F,row.names=F)
    write.table(ASEresult_1s_sum$geneOutput[f,],paste(infile,".ASE.1s.Gene.filter",sep=""),sep="\t",quote=F,row.names=F)
    write.table(ASEresult_1s_sum$geneOutput[f2,],paste(infile,".ASE.1s.Gene.filter.isoform",sep=""),sep="\t",quote=F,row.names=F)

}else if(args$type==2)
{
    mySEData <- SummarizedExperiment(
      assays=list(
	    lociAllele1Counts=matrix(
		    c(
			    inTable[,"lociAllele1Counts_tumor"],
			    inTable[,"lociAllele1Counts_normal"]
		    ),
		    ncol=2,
		    dimnames=list(
			    names(inTable_rowData), 
			    c('tumor', 'normal')
		    )
	    ),
	    lociAllele2Counts=matrix(
		    c(
			    inTable[,"lociAllele2Counts_tumor"],
			    inTable[,"lociAllele2Counts_normal"]
		    ),
		    ncol=2,
		    dimnames=list(
			    names(inTable_rowData), 
			    c('tumor', 'normal')
		    )
	    )
      ),
      rowData=inTable_rowData
    )
    ASEresults <- runMBASED(
    	ASESummarizedExperiment=mySEData,
    	isPhased=FALSE,
    	numSim=10^6,
    	BPPARAM = MulticoreParam()
    )
    ASEresult_2s_sum <- summarizeASEResults_2s(ASEresults)
    f  <- with(ASEresult_2s_sum$geneOutput, abs(majorAlleleFrequencyDifference) > 0.2 & pValueASE.adj < 0.05)
    f2 <- with(ASEresult_2s_sum$geneOutput, abs(majorAlleleFrequencyDifference) > 0.2 & pValueASE.adj < 0.05 & !is.na(pValueHeterogeneity) & pValueHeterogeneity <0.01 )
    write.table(ASEresult_2s_sum$geneOutput,paste(infile,".ASE.2s.Gene",sep=""),sep="\t",quote=F,row.names=F)
    write.table(ASEresult_2s_sum$geneOutput[f,],paste(infile,".ASE.2s.Gene.filter",sep=""),sep="\t",quote=F,row.names=F)
    write.table(ASEresult_2s_sum$geneOutput[f2,],paste(infile,".ASE.2s.Gene.filter.isoform",sep=""),sep="\t",quote=F,row.names=F)
}else
{
    stop("Invalidate -t/--type parameters. You must specify 1 OR 2")
}

save.image(paste(infile,".RData",sep=""))
