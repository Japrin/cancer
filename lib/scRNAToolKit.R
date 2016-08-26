loginfo <- function(msg) {
  timestamp <- sprintf("%s", Sys.time())
  msg <- paste0("[",timestamp, "] ", msg,"\n")
  cat(msg)
}

entrezToXXX<-function(x,type="SYMBOL",species="human")
{
  ret=c()
  x=as.character(x)
  if(species=="human")
  {
      suppressPackageStartupMessages(require("org.Hs.eg.db"))
  }else if(species=="mouse")
  {
      suppressPackageStartupMessages(require("org.Mm.eg.db"))
  }else
  {
      stop("species must be human or mouse\n")
  }
  if(species=="human")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Hs.egSYMBOL))[x]
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Hs.egENSEMBL))[x]
      }else if(type=="GENENAME")
      {
        ret=unlist(as.list(org.Hs.egGENENAME))[x]
      }
  }else if(species=="mouse")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Mm.egSYMBOL))[x]
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Mm.egENSEMBL))[x]
      }
  }
  return(as.character(ret))
}

XXXToEntrez<-function(x,type="SYMBOL",species="human")
{
  ret=c()
  x=as.character(x)
  if(species=="human")
  {
      suppressPackageStartupMessages(require("org.Hs.eg.db"))
  }else if(species=="mouse")
  {
      suppressPackageStartupMessages(require("org.Mm.eg.db"))
  }else
  {
      stop("species must be human or mouse\n")
  }
  if(species=="human")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Hs.egSYMBOL2EG))[x]  
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Hs.egENSEMBL2EG))[x]
      }
  }else if(species=="mouse")
  {
      if(type=="SYMBOL")
      {
        ret=unlist(as.list(org.Mm.egSYMBOL2EG))[x]  
      }else if(type=="ENSG")
      {
        ret=unlist(as.list(org.Mm.egENSEMBL2EG))[x]
      }
  }
  return(as.character(ret))
}

run.SIMLR <- function(inputData,n.clusters=4,myseed=NULL,my.title="",out.prefix,legend.txt,col.points,col.legend,
                      pch=16,width.pdf=10,height.pdf=8,margin.r=10,legend.inset=-0.25,legend.txt2=NULL,...){
    require(Matrix)
    require(parallel)
    require(igraph)
    # load the palettes for the plots
    require(grDevices)
    # load the SIMLR R package
    SIMLR.DIR <- "/Share/BP/zhenglt/01.bin/scRNASeq/SIMLR"
    source(sprintf("%s/R/SIMLR.R",SIMLR.DIR))
    source(sprintf("%s/R/compute.multiple.kernel.R",SIMLR.DIR))
    source(sprintf("%s/R/network.diffusion.R",SIMLR.DIR))
    source(sprintf("%s/R/utils.simlr.R",SIMLR.DIR))
    source(sprintf("%s/R/tsne.R",SIMLR.DIR))
    # load the C file
    dyn.load(sprintf("%s/src/projsplx_R.so",SIMLR.DIR))
    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%d) for SIMLR\n",myseed))

    res.SIMLR <- SIMLR(X=inputData,c=n.clusters)

    pdf(file=sprintf("%s.SIMLR.default.pdf",out.prefix),width=width.pdf,height=height.pdf)
    par(mar=c(5,5,4,margin.r),cex.lab=1.5,cex.main=1.5)
    ### color by predefined sample type
    plot(res.SIMLR$ydata, t='n', main=my.title,xlab="SIMLR component 1",ylab="SIMLR component 2")
    points(res.SIMLR$ydata,col=col.points,pch=pch,...)
    legend("right",legend=legend.txt,fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=16,border =NA,col = col.legend)
    if(!is.null(legend.txt2)) { 
           legend("topright",legend=legend.txt2,fill = NULL,inset = c(legend.inset,0),xpd = NA,cex=1.5,pch=unique(pch),border =NA,col = col.legend[1]) 
    }

    ## color by clustering result
    clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(res.SIMLR$y$cluster))),
                              names=as.character(unique(sort(res.SIMLR$y$cluster))))
    ##clusterColor["0"] <- "gray"
    print(clusterColor)
    plot(res.SIMLR$ydata, col=clusterColor[as.character(res.SIMLR$y$cluster)], 
         pch=16,main=my.title,xlab="SIMLR component 1",ylab="SIMLR component 2")
    legend("right",legend=sprintf("cluster%s",names(clusterColor)),
           fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=16,border =NA, 
           col = clusterColor)
    dev.off()


    return(res.SIMLR)
}


run.limma.from.matrixFile <- function(infile,designFile,fdr=0.1)
{
	inTable <- read.table(infile,row.names=1,check.names = F,header = T,stringsAsFactors = F)
	inMatrix <- as.matrix(inTable[,-1])

	suppressPackageStartupMessages(require(limma))
	suppressPackageStartupMessages(require("BiocParallel"))

	myDesign<-read.table(designFile,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
	design <- model.matrix(~patient+sampleType,data=myDesign)
	colnames(design)[length(colnames(design))]<-"II"

	register(MulticoreParam(4))
	fit <- lmFit(inMatrix, design)
	fit <- eBayes(fit)
	all.table  <- topTable(fit, coef = "II", n = Inf, sort = "p", p = 1)
	out.resSig <- topTable(fit, coef = "II", n = Inf, sort = "p", p = fdr)
	
	all.table <- cbind(data.frame(Gene=rownames(all.table)),
	      data.frame(Entrez=XXXToEntrez(rownames(all.table))),
	      all.table)
	
	out.resSig <- cbind(data.frame(Gene=rownames(out.resSig),stringsAsFactors=F),
	      data.frame(Entrez=XXXToEntrez(rownames(out.resSig)),stringsAsFactors=F),
	      out.resSig)

	#write.table(all.table,paste0(output.prefix,".all"),row.names = F,quote = F,sep = "\t")
	#write.table(out.resSig,paste0(output.prefix,".sig"),row.names = F,quote = F,sep = "\t")
	list(all=all.table,sig=out.resSig)
	#head(all.table)
	#head(out.resSig)
}

readGMT<-function(file) 
{
    f <- readLines(file)
    lst = sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
    names(lst) = sapply(lst, function(x) x[1])
    gSet = lapply(lst, function(x) x[-(1:2)])
    gLink = unlist(lapply(lst, function(x) x[2]))
    list(gSet=gSet,gLink=gLink)
}

gageAnalysis<-function(exp.fc,outDir,fdr=0.1)
{
    suppressPackageStartupMessages(require("gage"))
    suppressPackageStartupMessages(require("pathview"))
    suppressPackageStartupMessages(require("ReportingTools"))
    suppressPackageStartupMessages(require("hwriter"))
    data(kegg.gs)
    #c2.cp.gs=readGMT("/lustre1/zeminz_pkuhpc/00.database/MSigDB/msigdb_v4.0_files_to_download_locally/msigdb_v4.0_GMTs/c2.cp.v4.0.entrez.gmt")
    c2.cp.gs=readGMT("/Share/BP/zhenglt/00.database/MSigDB/msigdb_v4.0_files_to_download_locally/msigdb_v4.0_GMTs/c2.cp.v4.0.entrez.gmt")
    #### gage using kegg.gs
    gagePathway<-function(exp.fc,outDir,gs,gs.name,fdr,gl=NULL)
    {
	# 1 direction
	fc.p <- gage(exp.fc, gsets = gs, ref = NULL, samp = NULL)
	fc.p$greater=as.data.frame(fc.p$greater)
	fc.p$less=as.data.frame(fc.p$less)
	sel.h <- fc.p$greater[, "q.val"] < fdr & !is.na(fc.p$greater[, "q.val"])
	path.ids.h <- rownames(fc.p$greater)[sel.h]
	sel.l <- fc.p$less[, "q.val"] < fdr & !is.na(fc.p$less[,"q.val"])
	path.ids.l <- rownames(fc.p$less)[sel.l]
	# 2 direction
	fc.p.2d <- gage(exp.fc, gsets = gs, ref = NULL, samp = NULL, same.dir = F)
	fc.p.2d$greater=as.data.frame(fc.p.2d$greater)
	fc.p.2d$less=as.data.frame(fc.p.2d$less)
	sel.2d <- fc.p.2d$greater[, "q.val"] < fdr & !is.na(fc.p.2d$greater[, "q.val"])
	path.ids.2d <- rownames(fc.p.2d$greater)[sel.2d]
	# all path ids
	path.ids <- c(path.ids.h, path.ids.l, path.ids.2d)
	if(sum(regexpr("hsa[0-9]{5}",names(gs),perl=T)==rep(1,length(names(gs))))==length(names(gs)))
	{
		path.ids <- substr(path.ids, 1, 8)
	}
	## output table
	out1<-fc.p$greater[sel.h,]
	out1<-data.frame(pathway=rownames(out1),out1)
	out2<-fc.p$less[sel.l,]
	out2<-data.frame(pathway=rownames(out2),out2)
	out1d<-rbind(out1,out2)
	out2d<-fc.p.2d$greater[sel.2d,]
	out2d<-data.frame(pathway=rownames(out2d),out2d)
	write.table(out1d, file=paste(outDir,"/gage.",gs.name,".1d.txt",sep=""),sep="\t",quote=F,row.names=F)
	write.table(out2d, file=paste(outDir,"/gage.",gs.name,".2d.txt",sep=""),sep="\t",quote=F,row.names=F)
	####
	if(gs.name=="KEGG")
	{
	    # pathview, only for "KEGG"
	    picDir <- paste(outDir,"/reports/figure_gage",sep="")
	    dir.create(picDir,recursive=T, showWarnings = FALSE)
	    out.suffix="gage"
	    oriDir <- getwd()
	    setwd(picDir)
	    ## OUT/reports/figure_gage/hsa04740.gage.DESeq2.png
	    #pv.out.list <- sapply(path.ids, function(pid) pathview(kegg.dir="/lustre1/zeminz_pkuhpc/00.database/kegg/pathview",gene.data =  exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))
	    pv.out.list <- sapply(path.ids, function(pid) pathview(kegg.dir="/Share/BP/zhenglt/00.database/kegg/pathview",gene.data =  exp.fc, pathway.id = pid, species = "hsa", out.suffix=out.suffix))
	    setwd(oriDir)
	    if(nrow(out1d) > 0)
	    {
		imagename <- paste0("figure_gage/",substr(out1d$pathway, 1, 8),".",out.suffix,".png")
		out1d$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=50, width=50)
	    }
	    if(nrow(out2d) > 0)
	    {
		imagename <- paste0("figure_gage/",substr(out2d$pathway, 1, 8),".",out.suffix,".png")
		out2d$Image <- hwriteImage(imagename, link = imagename, table = FALSE, height=50, width=50)
	    }
	}
	#### Report writting
	addGSetLink <- function(object, ...)
	{
		if(!is.null(gl))
		{
			object$pathway <- hwrite(as.character(object$pathway), link=gl[as.character(object$pathway)], table = FALSE)
		}
		return(object)
	}
	gageReport <- HTMLReport(shortName = paste("gage_analysis_rnaseq_",gs.name,"",sep=""), title = paste("GAGE analysis of RnaSeqData (", gs.name, ")",sep=""), reportDirectory = "reports", basePath=outDir)
	publish(hwrite("One direction (up or down)", heading=2), gageReport)
	if(nrow(out1d) > 0)
	{
		publish(out1d, gageReport, reportDir="./reports", .modifyDF=list(addGSetLink))
	}
	publish(hwrite("Two direction (up and down)", heading=2), gageReport)
	if(nrow(out2d) > 0)
	{
		publish(out2d, gageReport, reportDir="./reports", .modifyDF=list(addGSetLink))
	}
	finish(gageReport)
	list(report=gageReport)
    }
    gage.kegg=gagePathway(exp.fc,outDir,kegg.gs,"KEGG",0.1)
    gage.Canonical=gagePathway(exp.fc,outDir,c2.cp.gs$gSet,"Canonical",fdr,c2.cp.gs$gLink)
    ## return report obj for "index" web
    list(gage.kegg$report, gage.Canonical$report)
}

outputDEGene<-function(res,resSigStrict,dds,vstMat,fc.df,outDir)
{
    ## txt report
    resSig <- subset(res, padj < 0.1)
    resSigStrict <- subset(resSig, padj < opt.qval & abs(log2FoldChange) > 1)

    out_res <- as.data.frame(res)
    out_res$EntrezGeneID <- rownames(out_res)
    out_res$GeneSymbol <- entrezToXXX(out_res$EntrezGeneID)
    out_res$ENSG <- entrezToXXX(out_res$EntrezGeneID,"ENSG")
    
    if(nrow(fc.df)>0)
    {
    	out_res <- cbind(out_res, as.data.frame(vstMat)[rownames(out_res),], fc.df[rownames(out_res),c(-1,-2,-3)])
    }else
    {
    	out_res <- cbind(out_res, as.data.frame(vstMat)[rownames(out_res),])
    }
    
    out_resSig <- subset(out_res, padj < 0.1 )
    out_resSigStrict <- subset(out_resSig, padj < opt.qval & abs(log2FoldChange) > 1 )
    write.table(out_res, file=paste(outDir,"/DESeq2.txt",sep=""),sep="\t",quote=F,row.names=F)
    write.table(out_resSig, file=paste(outDir,"/DESeq2.sig.txt",sep=""),sep="\t",quote=F,row.names=F)
    write.table(out_resSigStrict, file=paste(outDir,"/DESeq2.sig.strict.txt",sep=""),sep="\t",quote=F,row.names=F)
    ## HTML report
    des2Report <- HTMLReport(shortName = 'RNAseq_analysis_with_DESeq2', title = 'RNA-seq analysis of differential expression using DESeq2', reportDirectory="reports", basePath=outDir)
    if(nrow(resSigStrict)>0)
    {
	    publish(resSigStrict,des2Report, pvalueCutoff=opt.qval,lfc=1, n=Inf, annotation.db="org.Hs.eg.db", DataSet=dds, factor = colData(dds)$sampleType, reportDir="./reports")
    }
    finish(des2Report)
    des2Report
}


## resSig: result return by run.limma.from.matrixFile
outputDEGeneFromLimma <-function(res,output.prefix)
{
    ## txt report
    all.table <- res$all
    resSig <- res$sig
    #resSigStrict <- subset(resSig, adj.P.Val < 0.05 & abs(logFC) > 1 )
    write.table(all.table,paste0(output.prefix,".all.txt"),row.names = F,quote = F,sep = "\t")
    write.table(resSig,paste0(output.prefix,".sig.txt"),row.names = F,quote = F,sep = "\t")
    #write.table(resSigStrict,paste0(output.prefix,".sig.strict.txt"),row.names = F,quote = F,sep = "\t")

    suppressPackageStartupMessages(require("ReportingTools"))
    suppressPackageStartupMessages(require("hwriter"))
    ## HTML report
    deReport <- HTMLReport(shortName = 'RNA_analysis_with_Limma', title = 'RNA Analysis of differential expression using Limma', reportDirectory="reports", basePath=dirname(output.prefix))
    addEntrezLink <- function(object, ...)
    {
    	object$Entrez <- hwrite(as.character(object$Entrez), link = paste0("http://www.ncbi.nlm.nih.gov/gene/", as.character(object$Entrez)), table = FALSE)
	return(object)
    }
    if(nrow(resSig)>0)
    {
    	publish(resSig,deReport, annotation.db="org.Hs.eg.db",.modifyDF=list(addEntrezLink) )
    }
    finish(deReport)
    deReport
}

hyperGEA<-function(x,selectedIDs,universeIDs,reportName,outDir, ...)
{
    suppressPackageStartupMessages(require("ReportingTools"))
    suppressPackageStartupMessages(require("hwriter"))
    suppressPackageStartupMessages(require("GOstats"))
    suppressPackageStartupMessages(require("Category"))
    suppressPackageStartupMessages(require("KEGG.db"))
    
    addGOIDLink <- function(object, ...)
    {
    	object$GOID <- hwrite(as.character(object$GOID), link = paste0("http://amigo.geneontology.org/amigo/term/", as.character(object$GOID)), table = FALSE)
	return(object)
    }
    addKEGGIDLink <- function(object, ...)
    {
    	#http://www.genome.jp/kegg/pathway/hsa/hsa05219.html
    	object$KEGGID <- hwrite(as.character(object$KEGGID), link = paste0("http://www.genome.jp/kegg/pathway/hsa/hsa", as.character(object$KEGGID),".html"), table = FALSE)
	return(object)
    }

    aReport <- HTMLReport(shortName = paste("Hypergeometric Tests for Gene Set (",reportName,")",sep=""), title = paste("Hypergeometric Tests for Gene Set (",reportName,")",sep=""), reportDirectory="reports", basePath=outDir)

    if(length(selectedIDs)>0)
    {
	aParams <- new(x, geneIds = selectedIDs, universeGeneIds = universeIDs, annotation ="org.Hs.eg", pvalueCutoff = 0.01, testDirection = "over", ... )
	aResults <- hyperGTest(aParams)
	aResultsSummary <- summary(aResults)

	if(nrow(aResultsSummary)>0)
	{
		if(x=="GOHyperGParams")
		{
		    #http://amigo.geneontology.org/amigo/term/GO:0000070
		    names(aResultsSummary)[1]<-"GOID"
		    publish(aResultsSummary, aReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg", pvalueCutoff= 0.01, .modifyDF=list(addGOIDLink))
		}else if(x=="KEGGHyperGParams")
		{
		    publish(aResultsSummary, aReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg", pvalueCutoff= 0.01, .modifyDF=list(addKEGGIDLink))
		}
	}
    }
    finish(aReport)
    aReport
   
    # not work still
    #pfamParams <- new("PFAMHyperGParams", geneIds= selectedIDs, universeGeneIds=universeIDs, annotation="org.Hs.eg", pvalueCutoff= 0.01, testDirection="over")
    #PFAMResults <- hyperGTest(pfamParams)
    #PFAMReport <- HTMLReport(shortName = 'pfam_analysis_rnaseq', title = "PFAM analysis of RnaSeqData", reportDirectory = paste(outDir,"/reports",sep=""))
    #publish(PFAMResults, PFAMReport, selectedIDs=selectedIDs, annotation.db="org.Hs.eg",categorySize=5)
    #finish(PFAMReport)

}

pairs2 <- function (x, labels, panel = points, ..., lower.panel = panel, 
    upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
    label.pos = 0.5 + has.diag/3, line.main = 3, cex.labels = NULL, 
    font.labels = 1, row1attop = TRUE, gap = 1, log = "",xlim=NULL, ylim=NULL) 
{
    if (doText <- missing(text.panel) || is.function(text.panel)) 
        textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
            y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
        oma, ...) {
        xpd <- NA
        if (side%%2L == 1L && xl[j]) 
            xpd <- FALSE
        if (side%%2L == 0L && yl[i]) 
            xpd <- FALSE
        if (side%%2L == 1L) 
            Axis(x, side = side, xpd = xpd, ...)
        else Axis(y, side = side, xpd = xpd, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]])) 
                x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]]))) 
                stop("non-numeric argument to 'pairs'")
        }
    }
    else if (!is.numeric(x)) 
        stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
        stop("only one column in the argument to 'pairs'")
    if (doText) {
        if (missing(labels)) {
            labels <- colnames(x)
            if (is.null(labels)) 
                labels <- paste("var", 1L:nc)
        }
        else if (is.null(labels)) 
            doText <- FALSE
    }
    oma <- if ("oma" %in% nmdots) 
        dots$oma
    main <- if ("main" %in% nmdots) 
        dots$main
    if (is.null(oma)) 
        oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    xl <- yl <- logical(nc)
    if (is.numeric(log)) 
        xl[log] <- yl[log] <- TRUE
    else {
        xl[] <- grepl("x", log)
        yl[] <- grepl("y", log)
    }
    for (i in if (row1attop) 
        1L:nc
    else nc:1L) for (j in 1L:nc) {
        l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", 
            ""))
        # modify here
         if (is.null(xlim) & is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l)
	if (is.null(xlim) & !is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l, ylim=ylim)
	if (!is.null(xlim) & is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l, xlim = xlim)
	if (!is.null(xlim) & !is.null(ylim))
		localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
			type = "n", ..., log = l, xlim = xlim, ylim=ylim)
        #localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
        #    type = "n", ..., log = l)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2L) || !has.upper || !has.lower)) 
                localAxis(1L + 2L * row1attop, x[, j], x[, i],
                  ...)
            # modify here
	    #if (i == nc && (j%%2L || !has.upper || !has.lower)) 
            #    localAxis(3L - 2L * row1attop, x[, j], x[, i], 
            #      ...)
            #if (j == 1 && (!(i%%2L) || !has.upper || !has.lower)) 
            #    localAxis(2L, x[, j], x[, i], ...)
            if (j == nc && (i%%2L || !has.upper || !has.lower)) 
                localAxis(4L, x[, j], x[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  localDiagPanel(as.vector(x[, i]), ...)
                if (doText) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  xlp <- if (xl[i]) 
                    10^0.5
                  else 0.5
                  ylp <- if (yl[j]) 
                    10^label.pos
                  else label.pos
                  text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                    font = font.labels)
                }
            }
            else if (i < j) 
                localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                  i]), ...)
            else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                i]), ...)
            if (any(par("mfg") != mfg)) 
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
            font = font.main)
    }
    invisible(NULL)
}

panel.cor <- function (x, y, digits = 2, meth = "pearson", cex.cor = 1,...)
{
	#print(meth)
	usr <- par("usr")
	on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
        r <- round(cor(x, y, method = meth), digits)
	cr <- round(cor.test(x, y, method = meth, alternative = "t")$p.value,digits)
	text(0.5, 0.65, paste("r =", r), cex = cex.cor)
	if (!cr == 0) {
		text(0.5, 0.35, paste("p =", cr), cex = cex.cor)
	}
	else text(0.5, 0.35, paste("p < 0.01"), cex = cex.cor)
}


#' group specific genes
#' 
#' @param dat.g
#' @param grps
#' @param out.prefix (string) output prefix
#' @param FDR.THRESHOLD
#' @param FC.THRESHOLD
runMultiGroupSpecificGeneTest <- function(dat.g,grps,out.prefix,mod=NULL,FDR.THRESHOLD=0.05,FC.THRESHOLD=1,verbose=F,n.cores=NULL)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("doParallel"))
    g <- unique(grps)
    dat.g <- as.matrix(dat.g)
    registerDoParallel(cores = n.cores)
    ret <- ldply(rownames(dat.g),function(v){
                  aov.out <- aov(y ~ g,data=data.frame(y=dat.g[v,],g=grps))
                  aov.out.s <- summary(aov.out)
                  t.res.f <- unlist(aov.out.s[[1]]["g",c("F value","Pr(>F)")])
                  aov.out.hsd <- TukeyHSD(aov.out)
                  hsd.name <- rownames(aov.out.hsd$g)
                  t.res.hsd <- c(aov.out.hsd$g[,"diff"],aov.out.hsd$g[,"p adj"])
                  t.res.hsd.minP <- min(aov.out.hsd$g[,"p adj"])
                  j <- which.min(aov.out.hsd$g[,"p adj"])
                  t.res.hsd.minPDiff <- if(is.na(t.res.hsd.minP)) NaN else aov.out.hsd$g[j,"diff"]
                  t.res.hsd.minPCmp <- if(is.na(t.res.hsd.minP)) NaN else rownames(aov.out.hsd$g)[j]
                  ## whether cluster specific ?
                  t.res.spe  <-  sapply(g,function(v){ all( aov.out.hsd$g[grepl(v,hsd.name,perl=T),"p adj"] < FDR.THRESHOLD )  })
                  ## wheter up across all comparison ?
                  is.up <- sapply(g,function(v){  all( aov.out.hsd$g[grepl(paste0(v,"-"),hsd.name),"diff"]>0 ) & all( aov.out.hsd$g[grepl(paste0("-",v),hsd.name),"diff"]<0 ) }) 
                  is.down <- sapply(g,function(v){  all( aov.out.hsd$g[grepl(paste0(v,"-"),hsd.name),"diff"]<0 ) & all( aov.out.hsd$g[grepl(paste0("-",v),hsd.name),"diff"]>0 ) })
                  is.clusterSpecific <- (sum(t.res.spe,na.rm = T) == 1)
                  if(is.clusterSpecific){
                    t.res.spe.lable <- names(which(t.res.spe))  
                    if(is.up[t.res.spe.lable]) { 
                        t.res.spe.direction <- "UP" 
                    }else if(is.down[t.res.spe.lable]) {
                        t.res.spe.direction <- "DOWN"
                    }else{
                        t.res.spe.direction <- "INCONSISTANT"
                    }
                  }else{
                    t.res.spe.lable <- "NA"
                    t.res.spe.direction <- "NA"
                  }
                  if(!is.null(mod) && mod=="cluster.specific") { }
                  structure(c(t.res.f,t.res.hsd,t.res.hsd.minP,t.res.hsd.minPDiff,t.res.hsd.minPCmp,
                              t.res.spe,is.clusterSpecific,t.res.spe.lable,t.res.spe.direction),
                            names=c("F","F.pvalue",paste0("HSD.diff.",hsd.name),paste0("HSD.padj.",hsd.name),
                                    "HSD.padj.min","HSD.padj.min.diff","HSD.padj.min.cmp",
                                    paste0("cluster.specific.",g),"is.clusterSpecific","cluster.lable","cluster.direction"))
                  
            },.progress = "none",.parallel=T)
    #print(str(ret))
    ret.df <- data.frame(geneID=rownames(dat.g),geneSymbol=entrezToXXX(rownames(dat.g)),stringsAsFactors=F)
    ret.df <- cbind(ret.df,ret)
    rownames(ret.df) <- rownames(dat.g)
    ## type conversion
    i <- 3:(ncol(ret.df)-(4+length(g)))
    ret.df[i]<-lapply(ret.df[i],as.numeric)
    i <- (ncol(ret.df)-(length(g)+2)):(ncol(ret.df)-2)
    ret.df[i]<-lapply(ret.df[i],as.logical)
    ## adjust F test's p value
    ret.df$q <- 1
    ret.df.1 <- subset(ret.df,!is.na(F.pvalue))
    ret.df.1$q <- p.adjust(ret.df.1[,"F.pvalue"],method = "BH")
    ret.df.2 <- subset(ret.df,is.na(F.pvalue))
    ret.df <- rbind(ret.df.1,ret.df.2)
    ret.df <- ret.df[order(ret.df$q,ret.df$HSD.padj.min),]
    ### select
    ret.df.sig <- subset(ret.df,q<FDR.THRESHOLD & HSD.padj.min<FDR.THRESHOLD & abs(HSD.padj.min.diff)>=FC.THRESHOLD)
    ### output
    write.table(ret.df.sig,file = sprintf("%s.aov.sig.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    if(verbose){
        write.table(ret.df,file = sprintf("%s.aov.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    }
    #print(str(ret.df))
    return(list(aov.out=ret.df,aov.out.sig=ret.df.sig))
}

#' run t-test on given data(row for genes and column for samples)
#' 
#' @param dat.g1
#' @param dat.g2
#' @param out.prefix (string) output prefix
#' @param FDR.THRESHOLD
runTTest <- function(dat.g1,dat.g2,out.prefix,FDR.THRESHOLD=0.05,FC.THRESHOLD=1,verbose=F,n.cores=NULL)
{
    suppressPackageStartupMessages(require("plyr"))
    suppressPackageStartupMessages(require("doParallel"))
    registerDoParallel(cores = n.cores)
    dat.g1 <- as.matrix(dat.g1)
    dat.g2 <- as.matrix(dat.g2)
    ret <- ldply(rownames(dat.g2),function(v){
           x <- dat.g1[v,]
           y <- dat.g2[v,]
           x.mean <- mean(x)
           y.mean <- mean(y)
           t.res.fc  <- y.mean-x.mean
           t.res.stat <- NA
           t.res.p <- 1
           if(t.res.fc!=0){
               t.res <- t.test(x=x,y=y)
               t.res.p <- t.res$p.value
               t.res.stat <- t.res$statistic
           }
           structure(c(t.res.p,y.mean,x.mean,t.res.fc,t.res.stat),
                     names=c("p.value","y.mean","x.mean","fc","stat"))

    },.progress = "none",.parallel=T)

    rownames(ret) <- rownames(dat.g2)
    ret.df <- data.frame(geneID=rownames(ret),geneSymbol=entrezToXXX(rownames(ret)),stringsAsFactors=F)
    ret.df <- cbind(ret.df,ret)
    ret.df$p.adj <- 1
    ret.df.1 <- subset(ret.df,fc!=0)
    ret.df.2 <- subset(ret.df,fc==0)
    ret.df.1$p.adj <- p.adjust(ret.df.1[,"p.value"],method = "BH")
    ret.df <- rbind(ret.df.1,ret.df.2)
    ret.df <- ret.df[order(ret.df$p.adj,ret.df$p.value),]
    ret.df.sig <- subset(ret.df,p.adj<FDR.THRESHOLD & abs(fc)>=FC.THRESHOLD)
    write.table(ret.df.sig,file = sprintf("%s.ttest.sig.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    if(verbose){
        write.table(ret.df,file = sprintf("%s.ttest.txt",out.prefix),quote = F,row.names = F,col.names = T,sep = "\t")
    }
    return(list(ttest.out=ret.df,ttest.out.sig=ret.df.sig))
}

#' run t-SNE analysis on given data(row for genes and column for samples), transpose the data internally
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param legend
#' @param col.points 
#' @param col.legend
runTSNEAnalysis <- function(dat.plot,out.prefix,legend,col.points,col.legend,pch=16,pch.legend=16,inPDF=TRUE,eps=2.0,dims=2,k=NULL,do.dbscan=F,myseed=NULL,width.pdf=10,height.pdf=8,margin.r=8,legend.inset=-0.21,preSNE=NULL,...)
{
    suppressPackageStartupMessages(require("Rtsne"))
    suppressPackageStartupMessages(require("dbscan"))
    suppressPackageStartupMessages(require("RColorBrewer"))
    if(is.null(k)) { k=dims+1; }

    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    dat.plot <- t(dat.plot)
    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%d) for tsne\n",myseed))
    ret.list <- list()
    doit <- function(par.perplexity)
    {
        if(inPDF){
            pdf(file=sprintf("%s.perplexity%d.eps%4.2f.minPts%d.dims%d.pdf",out.prefix,par.perplexity,eps,k,dims),width=width.pdf,height=height.pdf)
            par(mar=c(5,5,4,margin.r),cex.lab=1.5,cex.main=1.5)
        }
        loginfo(sprintf("... begin Rtsne\n"))
        if(is.null(preSNE)) { 
            Rtsne.res <- Rtsne(dat.plot,perplexity=par.perplexity,dims=dims) 
        }else{
            Rtsne.res <- preSNE
        }
        loginfo(sprintf("... end Rtsne\n"))
        if(do.dbscan){
            loginfo(sprintf("... begin dbscan\n"))
            kNNdistplot(Rtsne.res$Y, k = k)
            abline(h = eps, lty = 2, lwd=1.5)
            dbscan.res <- dbscan(Rtsne.res$Y, eps = eps, minPts = k)
            loginfo(sprintf("... end dbscan\n"))
            ## cluster color
            clusterColor <- structure(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(dbscan.res$cluster))),
                              names=as.character(unique(dbscan.res$cluster)))
            clusterColor["0"] <- "gray"
            print(clusterColor)

            plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2],
                 col=clusterColor[as.character(dbscan.res$cluster)],
                 pch=pch,main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
            legend("right",legend=sprintf("cluster%d",unique(dbscan.res$cluster)),
                   fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=pch,border =NA,
                   col = clusterColor[as.character(unique(dbscan.res$cluster))])
            ret.list[[as.character(par.perplexity)]] <<- list(Rtsne.res=Rtsne.res,dbscan.res=dbscan.res,myseed=myseed)
        }else {
            ret.list[[as.character(par.perplexity)]] <<- list(Rtsne.res=Rtsne.res,myseed=myseed)
        }
        plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2], t='n', main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
        points(Rtsne.res$Y,col=col.points,pch=pch,...)
        legend("right",legend=legend,fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=pch.legend,border =NA,col = col.legend)
        
        plot(Rtsne.res$Y[,1],Rtsne.res$Y[,2], t='n', main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
        points(Rtsne.res$Y,col=col.points,pch=16,...)
        legend("right",legend=legend,fill = NULL,inset = legend.inset,xpd = NA,cex=1.5,pch=16,border =NA,col = col.legend)
        if(inPDF){
            dev.off()
        }
    }
    #sapply(seq(5,50,5), function(i) { tryCatch(doit(i), error = function(e) e, finally = loginfo(sprintf("runTSNEAnalysis() finally with perplexity %d \n",i)) ) } )
    sapply(seq(30,30,5), function(i) { tryCatch(doit(i), error = function(e) e, finally = loginfo(sprintf("runTSNEAnalysis() finally with perplexity %d \n",i)) ) } )
    return(ret.list)
}

#' run NMF analysis on given data(row for genes and column for samples)
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param ann.col  column annatation of dat.plot
#' @param ann.colors  a list specify how to mapping ann.col to colors
runNMFAnalysis <- function(dat.plot,out.prefix,ann.col,ann.colors)
{
    suppressPackageStartupMessages(require("NMF"))

    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    doNMF <- function(dat.plot)
    {
        nmf.res <- nmf(dat.plot, c(2:10,12,14,16), nrun = 100, .opt = "vp8", seed = 123456)
        pdf(file=sprintf("%s.pdf",out.prefix),width=12,height=10)
        print(plot(nmf.res))
        ###nmf.options(grid.patch=TRUE)
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         consensusmap(x[[i]],
                                      annCol = ann.col,
                                      annColors = ann.colors,
                                      fontsize=20) 
                    }, nmf.res$fit))
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         basismap(x[[i]],main = paste0("Metagenes (Rank=",i+1,")"),fontsize=20) 
                    }, nmf.res$fit))
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         coefmap(x[[i]],main = paste0("Metagene contributions in each sample (Rank=",i+1,")"),
                                 annCol = ann.col,
                                 annColors = ann.colors,
                                 fontsize=20)
                    }, nmf.res$fit))
        dev.off()
        save(nmf.res,file = sprintf("%s.RData",out.prefix))
    }
    tryCatch(doNMF(dat.plot), 
             error = function(e) e, 
             finally = loginfo(sprintf("runNMFAnalysis()'s finally")) )
}


runKMeansAnalysis <- function(in.data,out.prefix,sampleType,colSet,k,B=100,nfeatures=NULL)
{
    suppressPackageStartupMessages(require("cluster"))
    if(!is.null(nfeatures)) {
        dat.plot <- in.data[,seq_len(nfeatures)]
    }
    k.res <- kmeans(dat.plot,k,iter.max=1000,nstart=50)
    print(table(sampleType,k.res$cluster))
    gskmn <- clusGap(dat.plot, FUN = kmeans, nstart = 50, K.max = 30, B = B)
        
    pdf(file=sprintf("%s.clusGap.pdf",out.prefix),width=10,height=8)
    par(mar=c(5,5,4,8),cex.lab=1.5)
    plot(gskmn)
    dev.off()

}
#runKMeansAnalysis(pca.res$x,sprintf("%s/%s.het.PCA",out.dir,sample.id), myDesign$sampleType,sampleTypeColor[names(sampleTypeColor) %in% unique(as.character(myDesign$sampleType))], k=9,nfeatures=100)

runSC3Analysis <- function(in.data,out.prefix,sampleType,colSet,do.log.scale=FALSE,n.cores=4)
{
    suppressPackageStartupMessages(require("SC3"))
	suppressPackageStartupMessages(require("RColorBrewer"))
    ### hijack the original code
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/SC3/shiny-funcs.R")
    source("/Share/BP/zhenglt/02.pipeline/cancer/lib/SC3/iwanthue.R")

    cnames <- entrezToXXX(rownames(in.data))
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- rownames(in.data)[cnames.na]
    rownames(in.data) <- cnames

    ### helper functions
    output.cluster.labels <- function(values,input.param,kk)
    {
        clusts <- cutree(values$hc, kk)
        cell.order <- order.dendrogram(as.dendrogram(values$hc))
        d <- input.param$dataset
        colnames(d) <- clusts
        d <- d[ , cell.order]
        values$original.labels <- input.param$cell.names[cell.order]
        ## prepare a consensus matrix for downloading
        #tmp <- data.frame(values$consensus[cell.order, cell.order])
        #colnames(tmp) <- values$original.labels
        #rownames(tmp) <- seq(length=nrow(tmp))
        #values$consensus.download <- tmp

        values$new.labels <- reindex_clusters(colnames(d))
        colnames(d) <- values$new.labels
        values$dataset <- d
        values$col.gaps <- which(diff(as.numeric(colnames(d))) != 0)
        values$cell.labels <- data.frame(new.labels = as.numeric(values$new.labels), 
                                        original.labels = values$original.labels, 
                                        stringsAsFactors = FALSE)
        rownames(values$cell.labels) <- values$cell.labels$original.labels
        return(values)
    }

    plot.marker.genes <- function(out.prefix,values,input.param,input.threhold,kk,col.args.clusterAndCellType)
    {
        suppressPackageStartupMessages(require("ROCR"))
        d <- values$dataset
        col.gaps <- values$col.gaps

        mark.res <- get_marker_genes(d,
                                     as.numeric(colnames(d)),
                                     as.numeric(input.threhold$auroc.threshold),
                                     as.numeric(input.threhold$p.val.mark))
        if(is.null(mark.res)) {
            cat(sprintf("Unable to find significant marker genes from obtained clusters (k=%s).\n",kk))
            values$marker.genes <- data.frame(new.labels=c(),gene=c(),AUROC=c(),p.value=c(),stringsAsFactors = F)
            return(values)
        }
        ##colnames(mark.res) <- c("AUC","clusts","p.value")
        d.param <- mark_gene_heatmap_param(mark.res, unique(colnames(d)))
        #sc3.res$mark <- TRUE
        values$marker.genes <- data.frame(new.labels = as.numeric(mark.res[,2]),
                                        gene = rownames(mark.res),
                                        AUROC = as.numeric(mark.res[,1]),
                                        p.value = as.numeric(mark.res[,3]),
                                        stringsAsFactors = FALSE)
        write.table(values$marker.genes,sprintf("%s.markerGene.txt",out.prefix),sep = "\t",row.names = F,quote = F)

        dat.plot <- d[rownames(d.param$mark.res.plot), , drop = FALSE]
        colnames(dat.plot) <- as.character(values$hc$order)
        pheatmap::pheatmap(dat.plot,
                         show_colnames = FALSE,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         annotation_row = d.param$row.ann,
                         annotation_col = col.args.clusterAndCellType$col.ann.df, 
                         annotation_colors = col.args.clusterAndCellType$ann.colors,
                         annotation_names_row = FALSE,
                         gaps_row = d.param$row.gaps, gaps_col = col.gaps,
                         annotation_legend = FALSE,
                         fontsize = 20, fontsize_col = 12*55/ncol(d),fontsize_row = 250/max(nrow(d.param$mark.res.plot),25),
                         cellheight = 10)
        return(values)
    }

    output.outlier <- function(out.prefix,values,input.param,kk,sampleType)
    {
        suppressPackageStartupMessages(require("colorspace"))
        suppressPackageStartupMessages(require("ggplot2"))
        suppressPackageStartupMessages(require("rrcov"))
        # compute outlier cells
        #sc3.res$outl <- FALSE
        d <- values$dataset
        col.gaps <- values$col.gaps
        outl.res <- outl_cells_main(d, input.param$chisq.quantile)
        t <- as.data.frame(outl.res)
        colnames(t)[1] <- "outl"
        t$Cluster <- names(outl.res)
        t$Cells <- 1:dim(t)[1]
        t$Cluster <-  factor(t$Cluster, levels = unique(as.character(sort(as.numeric(t$Cluster)))))

        cols <- iwanthue(length(unique(t$Cluster)))
        Cells <- outl <- Cluster <- NULL

        #sc3.res$outl <- TRUE
        values$cells.outliers <- data.frame(new.labels = as.numeric(names(outl.res)),
                                          original.labels = values$original.labels,
                                          MCD.dist = outl.res,
                                          cell.type = sampleType[values$hc$order],
                                          stringsAsFactors = FALSE)
        rownames(values$cells.outliers) <- values$original.labels
        write.table(values$cells.outliers,sprintf("%s.outlier.txt",out.prefix),sep = "\t",row.names = F,quote = F)

        print(ggplot(t, aes(x = t$Cells, y = t$outl, fill = t$Cluster, color = t$Cluster)) +
            geom_bar(stat = "identity") +
            geom_point() +
            scale_fill_manual(values = cols) +
            scale_color_manual(values = cols) +
            guides(color = FALSE, fill = FALSE) +
            labs(y = "Outlier score") +
            theme_bw())
        ###ggsave(filename=sprintf("%s.outlier.pdf",out.prefix),width=8,height=6)
        return(values)
    }
    plot.DE.genes <- function(out.prefix,values,input.param,kk,sampleType,col.args.clusterAndCellType,FDR.THRESHOLD=0.05,FC.THRESHOLD=1,n.cores=n.cores)
    {
        dat.g <- values$dataset
        d.original.labels <- values$original.labels
        d.new.labels <- values$new.labels
        colnames(dat.g) <- d.original.labels
        sampleType <- sampleType[values$hc$order]
        names(sampleType) <- d.original.labels
        grps.list <- c()
        cmp.grp <- sprintf("Grp%02d",d.new.labels)
        for(j in unique(d.new.labels)){ 
            grps.list[[sprintf("Grp%02d",j)]] <- d.original.labels[d.new.labels==j]
        }
        annDF <- data.frame(cluster=structure(cmp.grp,names=colnames(dat.g)),stringsAsFactors = F)
        annDF[is.na(annDF)] <- "NA"
        annCol <- list(cluster=structure(auto.colSet(length(grps.list),name = "Dark2"),names=unique(cmp.grp)))
        if(length(unique(d.new.labels))>2){
            mgeneTest.out <- runMultiGroupSpecificGeneTest(dat.g,cmp.grp,
                                                           sprintf("%s.FDR%g",out.prefix,100*FDR.THRESHOLD),
                                                           FDR.THRESHOLD=FDR.THRESHOLD,
                                                           FC.THRESHOLD=FC.THRESHOLD, verbose=T,n.cores = n.cores)
            
            do.heatmap.plot <- function(g.list,extra="",kkk=1)
            {
                dat.g.sig <- dat.g[g.list,,drop=F]
                dat.tmp <- c()
                for(i in seq_along(grps.list))
                {
                    dat.tmp.1 <- dat.g.sig[,grps.list[[i]],drop=F]
                    if(ncol(dat.tmp.1)>=2){
                        dat.tmp <- cbind(dat.tmp,dat.tmp.1[,hclust(dist(t(dat.tmp.1)))$order,drop=F])
                    }else{
                        dat.tmp <- cbind(dat.tmp,dat.tmp.1)
                    }
                }
                dat.g.sig <- dat.tmp
                #annDF <- data.frame(cluster=structure(cmp.grp,names=colnames(dat.g.sig)),stringsAsFactors = F)
                #annDF[is.na(annDF)] <- "NA"
                #annCol <- list(cluster=structure(auto.colSet(length(grps.list),name = "Dark2"),names=unique(cmp.grp)))
                runHierarchicalClusteringAnalysis(dat.g.sig,mytitle = sprintf("K=%d",kk),
                                              pdf.width=18,pdf.height=15,do.clustering.col=F,
                                              sprintf("%s.aov.sig.noClusteringCol.FDR%g%s",
                                                      out.prefix,100*FDR.THRESHOLD,extra), 
                                              sampleType=sampleType[colnames(dat.g.sig)], 
                                              colSet=colSet,
                                              ann.extra.df = annDF,
                                              ann.extra.df.col = annCol,
                                              ann.bar.height = 0.8,
                                              k.row=kkk,clonotype.col=NULL,ntop=NULL, 
                                              row.names.original=TRUE,
                                              complexHeatmap.use=TRUE,verbose=FALSE)
            }
            ## aov genes
            g.list <- as.character(rownames(mgeneTest.out$aov.out.sig))
            if(length(g.list)>=3){ 
                do.heatmap.plot(g.list,extra="") 
                if(length(g.list)>30){ do.heatmap.plot(head(g.list,n=30),extra=".top30") }
            }
            ## cluster specific genes
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" )))
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific", kk=1) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" & cluster.lable==t.grp ))),n=10))
            }
            ###cat(sprintf("#### TEST kk=%d ####",kk))
            ###print(g.list)
            ###print(str(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction!="INCONSISTANT" & cluster.lable==cmp.grp[1])))
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.top10", kk=1) }
            ## cluster specific genes & up-regulated
            g.list <- as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="UP")))
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.UP", kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }
            g.list <- c()
            for(t.grp in unique(cmp.grp)){
                g.list <- append(g.list, head(as.character(rownames(subset(mgeneTest.out$aov.out.sig,is.clusterSpecific==TRUE & cluster.direction=="UP" & cluster.lable==t.grp))),n=10))
            }
            if(length(g.list)>=3){ do.heatmap.plot(g.list,extra=".cluster.specific.UP.top10", kk=length(unique(subset(mgeneTest.out$aov.out.sig[g.list,],cluster.lable!="NA")$cluster.lable))) }

        }else if(length(unique(d.new.labels))==2){
            if(length(grps.list[[1]])<3 || length(grps.list[[2]])<3)
            {
                loginfo(sprintf("Two few samples: ncol(dat.g1)=%d, ncol(dat.g2)=%d",length(grps.list[[1]]),length(grps.list[[1]])))
            }else{
                ttest.out <<- runTTest(dat.g[,grps.list[[1]]],dat.g[,grps.list[[2]]],
                                       sprintf("%s.%s.FDR%g",out.prefix,"ttest",100*FDR.THRESHOLD),
                                       FDR.THRESHOLD,FC.THRESHOLD,verbose=T,n.cores = n.cores)
                do.heatmap.plot <- function(g.list,extra="",kkk=1){
                    if(length(g.list)>=3){
                            dat.g.1 <- dat.g[g.list,grps.list[[1]],drop=F]
                            dat.g.2 <- dat.g[g.list,grps.list[[2]],drop=F]
                            dat.g.sig <- cbind(dat.g.1[,hclust(dist(t(dat.g.1)))$order,drop=F], 
                                           dat.g.2[,hclust(dist(t(dat.g.2)))$order,drop=F])
                            runHierarchicalClusteringAnalysis(dat.g.sig,mytitle = sprintf("K=%d",kk),
                                                              pdf.width=18,pdf.height=15,
                                                              sprintf("%s.%s.sig.FDR%g.noClusteringCol%s",
                                                                      out.prefix,"ttest",100*FDR.THRESHOLD,extra), 
                                                              do.clustering.col=F,
                                                              sampleType=sampleType[colnames(dat.g.sig)], 
                                                              colSet=colSet,
                                                              ann.extra.df = annDF,
                                                              ann.extra.df.col = annCol,
                                                              ann.bar.height = 0.8,
                                                              k.row=length(unique(ttest.out$ttest.out.sig$fc>0)),
                                                              clonotype.col=NULL,ntop=NULL,
                                                              row.names.original=FALSE,
                                                              complexHeatmap.use=TRUE,verbose=FALSE)
                            dat.g.sig.df <- data.frame(geneSymbol=rownames(dat.g.sig))
                            dat.g.sig.df <- cbind(dat.g.sig.df,dat.g.sig)
                            write.table(dat.g.sig.df, 
                                        file = sprintf("%s.%s.sig.FDR%g.noClusteringCol%s.txt",
                                                       out.prefix,"ttest",100*FDR.THRESHOLD,extra), 
                                        quote = F,sep = "\t",row.names = F,col.names = T)
                    }
                } 
                ## for genes diff significantly
                g.list.1 <- as.character(rownames(ttest.out$ttest.out.sig))
                do.heatmap.plot(g.list.1,extra="")
                if(length(g.list.1)>30){ do.heatmap.plot(head(g.list.1,n=30),extra=".top30") }
                ## 
                g.list.1 <- c()
                g.list.1 <- append(g.list.1,head(as.character(rownames(subset(ttest.out$ttest.out.sig,x.mean > y.mean))),n=10))
                g.list.1 <- append(g.list.1,head(as.character(rownames(subset(ttest.out$ttest.out.sig,y.mean > x.mean))),n=10))
                if(length(g.list.1)>3){ do.heatmap.plot(g.list.1,extra=".ClusterTop10") }
            }
        }
    }
    ### old version(SC3_0.99.25)
    sc3(in.data, ks = 2:9, cell.filter = FALSE,gene.filter = FALSE,log.scale = do.log.scale,interactivity=FALSE,svm.num.cells=25000,use.max.cores=FALSE,n.cores=n.cores)
    ### latest version
    ###sc3(in.data, ks = 2:10, cell.filter = FALSE,gene.filter = FALSE,log.scale = do.log.scale,interactivity=FALSE,svm.num.cells=2500,n.cores=n.cores)
    input.threhold=list(auroc.threshold=0.85,p.val.mark=0.05)
    ###save.image(file = sprintf("%s.tmp.RData",out.prefix))
    ret.list=list()

    #for(kk in 2:16)
    for(kk in 2:10)
    {
        cat(sprintf("kk==%s\n",kk))
        ## extreac result
        sc3.cons<-with(sc3.interactive.arg, cons.table[cons.table[,1]=="spearman" & cons.table[,2]=="spectral" &　cons.table[,3]==kk,4][[1]])
        sc3.res<-list(consensus=sc3.cons[[1]],labels=sc3.cons[[2]],hc=sc3.cons[[3]],silh=sc3.cons[[4]])

        ## output cluster lables  
        sc3.res<-output.cluster.labels(sc3.res,sc3.interactive.arg,kk)
        ## plot silhouette
        pdf(sprintf("%s.k%s.pdf",out.prefix,kk),width = 10,height = 8)
        plot(sc3.res$silh)
        ## plot consensus
        col.ann.df <- data.frame(CellType=sampleType,Cluster=as.character(sc3.res$cell.labels[colnames(in.data),"new.labels"]))
        rownames(col.ann.df) <- colnames(sc3.res$consensus)
        #col.ann.df <- col.ann.df[sc3.res$hc$order,,drop=F]
        #ann.colors <- list(CellType=colSet,Cluster=as.character(seq_len(kk)))
        if(kk<=12){
            ann.colors <- list(CellType=colSet,Cluster=structure(brewer.pal(max(kk,3),"Set3"), .Names=as.character(seq_len(max(kk,3)))))
        }else{
            ann.colors <- list(CellType=colSet,Cluster=structure(colorRampPalette(brewer.pal(12,"Set3"))(kk), .Names=as.character(seq_len(max(kk,3)))))
        }
        pheatmap::pheatmap(sc3.res$consensus,
                           cluster_rows = sc3.res$hc, cluster_cols = sc3.res$hc,
                           cutree_rows = kk, cutree_cols = kk,
                           annotation_col = col.ann.df, annotation_colors = ann.colors,
                           fontsize = 20, fontsize_col = 12*55/max(ncol(in.data),25),
                           labels_col = colnames(in.data),
                           drop_levels = TRUE,
                           show_rownames = FALSE, show_colnames = TRUE)
        ## plot marker genes
        sc3.res<-plot.marker.genes(sprintf("%s.k%s",out.prefix,kk),sc3.res,sc3.interactive.arg,input.threhold,kk,list(col.ann.df=col.ann.df,ann.colors=ann.colors))
        ## outerlier
        sc3.res<-output.outlier(sprintf("%s.k%s",out.prefix,kk),sc3.res,sc3.interactive.arg,kk,sampleType)
        
        dev.off()
        
        tryCatch(
            plot.DE.genes(sprintf("%s.k%s",out.prefix,kk),sc3.res,sc3.interactive.arg,kk,sampleType,
                      list(col.ann.df=col.ann.df,ann.colors=ann.colors),
                      FDR.THRESHOLD=0.05,FC.THRESHOLD=1,n.cores=8),  
            error = function(e) e, 
            finally = loginfo(sprintf("plot.DE() finally")) )
        ret.list[[kk]] <- sc3.res
    }
    return(ret.list)
}

#' run PCA analysis on given data(row for genes and column for samples)
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param sampleType  vector of sample type
#' @param colSet  named vector of colors, names is the sample type, values is the mapped colors
#' @param ntop  select the ntop genes with the largest variance
#' @param verbose
runPCAAnalysis <- function(dat.plot,out.prefix,sampleType,colSet,ntop=NULL,verbose=FALSE,do.clust=FALSE,do.scale=TRUE,myseed=NULL,...)
{
    #require("DESeq2")
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("factoextra"))
    suppressPackageStartupMessages(require("FactoMineR"))
    rowVar <- apply(dat.plot,1,var)
    #rowVar <- apply(vstMat,1,function(x){ var(x)/(mean(x)^2) } )
    
    dat.plot <- as.matrix(dat.plot[rowVar>0,,drop=F])
    print(dim(dat.plot))
    print(length(sampleType))
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    if(!is.null(ntop)) { 
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
        dat.plot <- dat.plot[select,]
        n <- nrow(dat.plot)
    }
    dat.plot.t <- t(dat.plot)
    cnames <- entrezToXXX(colnames(dat.plot.t))
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- colnames(dat.plot.t)[cnames.na]
    colnames(dat.plot.t) <- cnames
    if(n<3) { loginfo(paste0("Too few genes: n=",n)); return(NULL) }
    if(m<3) { loginfo(paste0("Too few samples: n=",m)); return(NULL)
    }
    #### PCA
    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%d) for PCA\n",myseed))

    #### prcomp
    ####pca <- prcomp(dat.plot.t,scale=do.scale)
    ####percentVar <- pca$sdev^2/sum(pca$sdev^2)
    ####pca.eig <- (pca$sdev)^2
    ####variance.in.percent <- pca.eig*100/sum(pca.eig)
    ####cumvar.in.percent <- cumsum(variance.in.percent)
    ####pca.eig.df <- data.frame(eig = pca.eig, variance = variance.in.percent, cumvariance = cumvar.in.percent)

    ### FactoMinR
    pca <- PCA(dat.plot.t,ncp=min(30,n,m),graph=F,scale.unit=T)
    loginfo(sprintf("PCA() done.\n"))
    res.desc <- dimdesc(pca, axes = c(1:10))
    loginfo(sprintf("dimdesc() done.\n"))
    invisible(sapply(names(res.desc),function(x){
                     oo <- data.frame(gene=rownames(res.desc[[x]][["quanti"]]))
                     oo <- cbind(oo,res.desc[[x]][["quanti"]])
                     write.table(oo,sprintf("%s.variable.association.%s",out.prefix,x),row.names = F,col.names = T,sep="\t",quote = F)
             }))
    pca.eig.df <- pca$eig

    pca.var <- get_pca_var(pca)
    pca.ind <- get_pca_ind(pca)

    d <- data.frame(geneSymbol=rownames(pca.var$contrib))
    d <- cbind(d,pca.var$contrib)
    ozfile <- gzfile(sprintf("%s.contribution.txt.gz",out.prefix),open = "w")
	write.table(d,file = ozfile,sep = "\t",col.names = T,row.names = F,quote = F)
    close(ozfile)
    
    d <- data.frame(sampleType=sampleType, sampleName = colnames(dat.plot))
    d <- cbind(d,pca$ind$coord)
    ozfile <- gzfile(sprintf("%s.txt.gz",out.prefix),open = "w")
	write.table(d,file = ozfile,sep = "\t",col.names = T,row.names = F,quote = F)
    close(ozfile)
    
    #print(head(percentVar)) 
    nDim.plot <- min(30,nrow(pca.eig.df))
    pdf(file=sprintf("%s.pdf",out.prefix),width=10,height=8)
    par(mar=c(5,5,4,8),cex.lab=1.5)
    ### scree plot
    xat <- barplot(pca.eig.df[1:nDim.plot, 2], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")
    abline(h = 1*100/sum(pca.eig.df[,2]), lty = 2, col = "red", lwd=2)
    barplot(pca.eig.df[1:nDim.plot, 3], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Cumulative Variances", col ="steelblue")
    barplot(pca.eig.df[1:nDim.plot, 1], names.arg=1:nDim.plot, main = "Eigenvalue", xlab = "Principal Components", ylab = "Eigenvalues", col ="steelblue")

    ### loading plot
    my.plot.loading <- function(ld.x,ld.y,ld.x.lab,ld.y.lab){
        a <- seq(0, 2*pi, length = 100)
        plot( cos(a), sin(a), type = 'l', col="gray", xlab = ld.x.lab,  ylab = ld.y.lab)
        abline(h = 0, v = 0, lty = 2)
        arrows(0, 0, ld.x, ld.y, length = 0.1, angle = 15, code = 2)
        text(ld.x,ld.y, labels=names(ld.x), col="#0000FF99", cex = 1, adj=1)
        ####arrows(0, 0, pca.var$coord[, 1], pca.var$coord[, 2], length = 0.1, angle = 15, code = 2)
        ####text(pca.var$coord[, 1], pca.var$coord[, 2], labels=rownames(pca.var$coord), cex = 1, adj=1)
    }
    #my.plot.loading(pca$rotation[,1],pca$rotation[,2],"PC1","PC2")
    #my.plot.loading(pca$rotation[,1],pca$rotation[,3],"PC1","PC3")
    #my.plot.loading(pca$rotation[,2],pca$rotation[,3],"PC2","PC3")
    ### score plot
    my.plot.score <- function(pc.x,pc.y,pc.x.lab,pc.y.lab,pc.x.pve,pc.y.pve){
        plot(x=NULL,y=NULL, xlim = range(pc.x), ylim = range(pc.y), type = "n", cex.main=1.5,
             xlab=sprintf("%s: %s%% variance ",pc.x.lab, pc.x.pve), 
             ylab=sprintf("%s: %s%% variance ",pc.y.lab, pc.y.pve)) 
        invisible(lapply(levels(d$sampleType),function(x){ points(subset(d,sampleType==x,select=pc.x.lab)[,1],subset(d,sampleType==x,select=pc.y.lab)[,1],col=colSet[x],pch=16,cex=1.2)  }))
        legend("right",legend=names(colSet),fill = NULL,inset = -0.19,xpd = NA,cex=1.5,pch=16,border =NA,col = colSet)
    }
    my.plot.score(d$Dim.1,d$Dim.2,"Dim.1","Dim.2",round(pca.eig.df[1,2],2),round(pca.eig.df[2,2],2))
    my.plot.score(d$Dim.1,d$Dim.3,"Dim.1","Dim.3",round(pca.eig.df[1,2],2),round(pca.eig.df[3,2],2))
    my.plot.score(d$Dim.2,d$Dim.3,"Dim.2","Dim.3",round(pca.eig.df[2,2],2),round(pca.eig.df[3,2],2))

    print(fviz_screeplot(pca, ncp=nDim.plot))
    print(fviz_pca_var(pca, axes=c(1,2), col.var="contrib",geom=c("point","text")) + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    print(fviz_contrib(pca, choice = "var", axes = 1, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 1:2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 1, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 1:2, top = nDim.plot))
    #print(fviz_pca_var(pca, col.var="contrib") + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    #### extra components
    print(fviz_pca_var(pca, axes=c(3,4), col.var="contrib",geom=c("point","text")) + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    print(fviz_pca_var(pca, axes=c(5,6), col.var="contrib",geom=c("point","text")) + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    print(fviz_pca_var(pca, axes=c(3,6), col.var="contrib",geom=c("point","text")) + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    print(fviz_pca_var(pca, axes=c(4,5), col.var="contrib",geom=c("point","text")) + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    print(fviz_pca_var(pca, axes=c(1,5), col.var="contrib",geom=c("point","text")) + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    ###print(fviz_pca_ind(pca, axes=c(3,6),col.ind=colSet[sampleType],label="none"))
    print(fviz_contrib(pca, choice = "var", axes = 3, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 4, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 5, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 6, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = c(3,4), top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = c(5,6), top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = c(3,6), top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = c(4,5), top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = c(1,5), top = nDim.plot))
    #print(fviz_pca_ind(pca, label = "none", col.ind = patientcolors ,habillage = group, addEllipses = T) + theme_minimal())
    tryCatch(print(fviz_pca_ind(pca, label = "none", col.ind = colSet[sampleType],addEllipses = T) + theme_minimal()), error = function(e) e, finally = loginfo(sprintf("runPCAAnalysis() finally")) )
    dev.off()
    return(pca)
}

runPCAAnalysis.bak <- function(dat.plot,out.prefix,sampleType,colSet,ntop=NULL,verbose=FALSE,do.clust=FALSE,do.scale=TRUE,myseed=NULL,...)
{
    #require("DESeq2")
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("factoextra"))
    rowVar <- apply(dat.plot,1,var)
    #rowVar <- apply(vstMat,1,function(x){ var(x)/(mean(x)^2) } )
    
    dat.plot <- as.matrix(dat.plot)
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    if(!is.null(ntop)) { 
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
        dat.plot <- dat.plot[select,]
        n <- nrow(dat.plot)
    }
    dat.plot.t <- t(dat.plot)
    cnames <- entrezToXXX(colnames(dat.plot.t))
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- colnames(dat.plot.t)[cnames.na]
    colnames(dat.plot.t) <- cnames
    if(n<3) { loginfo(paste0("Too few genes: n=",n)); return(NULL) }
    if(m<3) { loginfo(paste0("Too few samples: n=",m)); return(NULL)
    }
    #### PCA
    if(is.null(myseed)){ myseed <- as.integer(Sys.time()) }
    set.seed(myseed)
    loginfo(sprintf("set.seed(%d) for PCA\n",myseed))
    pca <- prcomp(dat.plot.t,scale=do.scale)
    percentVar <- pca$sdev^2/sum(pca$sdev^2)

    pca.eig <- (pca$sdev)^2
    variance.in.percent <- pca.eig*100/sum(pca.eig)
    cumvar.in.percent <- cumsum(variance.in.percent)
    pca.eig.df <- data.frame(eig = pca.eig, variance = variance.in.percent, cumvariance = cumvar.in.percent)
    #eig.val <- get_eigenvalue(pca)
    pca.var <- get_pca_var(pca)
    pca.ind <- get_pca_ind(pca)

    d <- data.frame(geneSymbol=rownames(pca.var$contrib))
    d <- cbind(d,pca.var$contrib)
    ozfile <- gzfile(sprintf("%s.contribution.txt.gz",out.prefix),open = "w")
	write.table(d,file = ozfile,sep = "\t",col.names = T,row.names = F,quote = F)
    close(ozfile)
    
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], sampleType=sampleType, sampleName = colnames(dat.plot))
    ozfile <- gzfile(sprintf("%s.txt.gz",out.prefix),open = "w")
	write.table(d,file = ozfile,sep = "\t",col.names = T,row.names = F,quote = F)
    close(ozfile)
    
    #print(head(percentVar)) 
    nDim.plot <- min(30,nrow(pca.eig.df))
    pdf(file=sprintf("%s.pdf",out.prefix),width=10,height=8)
    par(mar=c(5,5,4,8),cex.lab=1.5)
    ### scree plot
    xat <- barplot(pca.eig.df[1:nDim.plot, "variance"], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")
    #lines(x = xat, pca.eig.df[1:nDim.plot, "variance"], type="b", pch=19, col = "red")
    abline(h = 1*100/sum(pca.eig), lty = 2, col = "red", lwd=2)
    xat <- barplot(pca.eig.df[1:nDim.plot, "cumvariance"], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Cumulative Variances", col ="steelblue")

    ### loading plot
    my.plot.loading <- function(ld.x,ld.y,ld.x.lab,ld.y.lab){
        a <- seq(0, 2*pi, length = 100)
        plot( cos(a), sin(a), type = 'l', col="gray", xlab = ld.x.lab,  ylab = ld.y.lab)
        abline(h = 0, v = 0, lty = 2)
        arrows(0, 0, ld.x, ld.y, length = 0.1, angle = 15, code = 2)
        text(ld.x,ld.y, labels=names(ld.x), col="#0000FF99", cex = 1, adj=1)
        ####arrows(0, 0, pca.var$coord[, 1], pca.var$coord[, 2], length = 0.1, angle = 15, code = 2)
        ####text(pca.var$coord[, 1], pca.var$coord[, 2], labels=rownames(pca.var$coord), cex = 1, adj=1)
    }
    #my.plot.loading(pca$rotation[,1],pca$rotation[,2],"PC1","PC2")
    #my.plot.loading(pca$rotation[,1],pca$rotation[,3],"PC1","PC3")
    #my.plot.loading(pca$rotation[,2],pca$rotation[,3],"PC2","PC3")
    ### score plot
    my.plot.score <- function(pc.x,pc.y,pc.x.lab,pc.y.lab,pc.x.pve,pc.y.pve){
        plot(x=NULL,y=NULL, xlim = range(pc.x), ylim = range(pc.y), type = "n", cex.main=1.5,
             xlab=sprintf("%s: %s%% variance ",pc.x.lab, pc.x.pve), 
             ylab=sprintf("%s: %s%% variance ",pc.y.lab, pc.y.pve)) 
        invisible(lapply(levels(d$sampleType),function(x){ points(subset(d,sampleType==x,select="PC1")[,1],subset(d,sampleType==x,select="PC2")[,1],col=colSet[x],pch=16,cex=1.2)  }))
        legend("right",legend=names(colSet),fill = NULL,inset = -0.19,xpd = NA,cex=1.5,pch=16,border =NA,col = colSet)
    }
    my.plot.score(d$PC1,d$PC2,"PC1","PC2",round(percentVar[1] * 100,2),round(percentVar[2] * 100,2))
    my.plot.score(d$PC1,d$PC3,"PC1","PC3",round(percentVar[1] * 100,2),round(percentVar[3] * 100,2))
    my.plot.score(d$PC2,d$PC3,"PC2","PC3",round(percentVar[2] * 100,2),round(percentVar[3] * 100,2))

    #library("FactoMineR")
    #pca <- PCA(dat.plot, graph = FALSE)
    print(fviz_screeplot(pca, ncp=nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 1, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "var", axes = 1:2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 1, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 2, top = nDim.plot))
    print(fviz_contrib(pca, choice = "ind", axes = 1:2, top = nDim.plot))
    print(fviz_pca_var(pca, col.var="contrib") + scale_color_gradient2(low="white", mid="blue", high="red", midpoint=50) + theme_minimal())
    #print(fviz_pca_ind(pca, label = "none", col.ind = patientcolors ,habillage = group, addEllipses = T) + theme_minimal())
    tryCatch(print(fviz_pca_ind(pca, label = "none", col.ind = patientcolors ,habillage = group, addEllipses = T) + theme_minimal()), error = function(e) e, finally = loginfo(sprintf("runPCAAnalysis() finally")) )
    dev.off()

    return(pca)
}


#' run hierachical clustering analysis on given data(row for genes and column for samples)
#' 
#' @param data.plot (matrix) input data (row for genes and column for samples)
#' @param out.prefix (string) output prefix
#' @param sampleType  vector of sample type
#' @param colSet  named vector of colors, names is the sample type, values is the mapped colors
#' @param k   number of clusters of the samples
#' @param clonotype.col   list for clonotype 
#' @param ntop  select the ntop genes with the largest variance
#' @param complexHeatmap.use  whether to use complexHeatmap
#' @param verbose
runHierarchicalClusteringAnalysis <- function(dat.plot,out.prefix,sampleType=NULL,colSet=NULL,k.row=1,k.col=1,
                                              clonotype.col=NULL,patient.col.list=NULL,ntop=NULL,
                                              complexHeatmap.use=FALSE,verbose=FALSE,
                                              row.names.original=FALSE,
                                              pdf.width=16,pdf.height=15,ann.extra.df=NULL,ann.extra.df.col=NULL,ann.bar.height=1.5,
                                              do.scale=TRUE,z.lo=-3,z.hi=3,z.step=1,z.title="Exp", annotation_legend_param=list(),
                                              do.cuttree=F,
                                              do.clustering.row=T,do.clustering.col=T,clustering.distance="spearman",clustering.method="complete",mytitle="",save.obj=F,...)
{
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("ComplexHeatmap"))
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("gridBase"))
    suppressPackageStartupMessages(require("dendextend"))
	suppressPackageStartupMessages(require("RColorBrewer"))
    
    dat.plot <- as.matrix(dat.plot)
    n <- nrow(dat.plot)
    m <- ncol(dat.plot)
    if(n<3) { loginfo(sprintf("Too few genes: n=%s",n)); return(NULL) }
    if(m<3) { loginfo(sprintf("Too few samples: m=%s",m)); return(NULL) }
    if(!is.null(ntop)) { 
        rowVar <- apply(dat.plot,1,var)
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
        dat.plot <- dat.plot[select,]
        n <- nrow(dat.plot)
    }
    if(!row.names.original){
        cnames <- entrezToXXX(rownames(dat.plot))
        cnames.na <- which(is.na(cnames))
        cnames[cnames.na] <- rownames(dat.plot)[cnames.na]
        rownames(dat.plot) <- cnames
    }

    ###rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
    if(verbose)
    {
        dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=if(!row.names.original) entrezToXXX(rownames(dat.plot)) else rownames(dat.plot) )
        dat.plot.export <- cbind(dat.plot.export,dat.plot)
        write.table(dat.plot.export,sprintf("%s.dat.txt",out.prefix),sep="\t",col.names = T,row.names = F,quote = F)
    }
    if(is.null(complexHeatmap.use))
    {
        patientcolors <- colSet[as.character(sampleType)]
        pdf(sprintf("%s.pdf",out.prefix),width=12,height=12)
        ####heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20), main="Most variable genes",cexRow=min(1.8,55/n),cexCol=min(1.8,55/m),distfun=function(x){ as.dist(1-cor(t(x),method = "spearman")) } )
        heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20),cexRow=min(1.8,55/n),cexCol=min(1.8,55/m),...)
        legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5)
        dev.off()
    }else
    {
        require("ComplexHeatmap")
        require("circlize")
        require("gridBase")
        require("dendextend")
        require("moduleColor")
        require("dynamicTreeCut")
        
        dat.plot.unscale <- dat.plot
        ###print(dat.plot[1:4,1:8])
        #### scale by row
        if(do.scale)
        {
            rowM <- rowMeans(dat.plot, na.rm = T)
            rowSD <- apply(dat.plot, 1, sd, na.rm = T)
            dat.plot <- sweep(dat.plot, 1, rowM)
            dat.plot <- sweep(dat.plot, 1, rowSD, "/")
            dat.plot[dat.plot < -3] <- -3
            dat.plot[dat.plot > 3] <- 3
            ###print(dat.plot[1:4,1:8])
        }else{
            tmp.var <- pretty(abs(dat.plot),n=5)
            z.lo <- tmp.var[1]
            z.hi <- tmp.var[length(tmp.var)]
            z.step <- tmp.var[2]-tmp.var[1]
        }
        branch.row <- FALSE
        branch.col <- FALSE
        obj.hclust.col <- NULL
        obj.hclust.row <- NULL
        if(do.clustering.col){
            if(clustering.distance=="spearman" || clustering.distance=="pearson"){
                tryCatch({
                        obj.hclust.col <- hclust(as.dist(1-cor(dat.plot.unscale,method=clustering.distance)),method=clustering.method)
                        branch.col <- color_branches(as.dendrogram(obj.hclust.col),k=k.col)
                    },
                    error = function(e) { cat("using spearman or pearson as distance failed; will try to fall back to use euler distance ... \n"); e }
                    ###error = function(e) { tmp.dat <<- dat.plot.unscale; print(clustering.method); print(k.col); e }
                    )
            }
            if(is.logical(branch.col) && !branch.col){
                obj.hclust.col <- hclust(dist(t(dat.plot.unscale)),method=clustering.method)
                branch.col <- color_branches(as.dendrogram(obj.hclust.col),k=k.col)
            }
        }
        if(do.clustering.row){
            if(clustering.distance=="spearman" || clustering.distance=="pearson"){
                tryCatch({
                        obj.hclust.row <- hclust(as.dist(1-cor(t(dat.plot.unscale),method=clustering.distance)),method=clustering.method)
                        branch.row <- color_branches(as.dendrogram(obj.hclust.row),k=k.row)
                    },
                    error = function(e) { cat("using spearman or pearson as distance failed; will try to fall back to use euler distance ... \n"); e }
                    )
            }
            if(is.logical(branch.row) && !branch.row){
                obj.hclust.row <- hclust(dist(dat.plot.unscale),method=clustering.method)
                branch.row <- color_branches(as.dendrogram(obj.hclust.row),k=k.row)
            }
        }
        if(do.cuttree && !is.null(obj.hclust.col)){
            ### cut tree
            dend.cutree <- cutree(obj.hclust.col, k=2:10, order_clusters_as_data = T)
            colnames(dend.cutree) <- sprintf("K=%s",colnames(dend.cutree))
            ### output clusters
            dend.cutree.df <- data.frame(sampleID=rownames(dend.cutree))
            dend.cutree.df <- cbind(dend.cutree.df,dend.cutree)
            write.table(dend.cutree.df,sprintf("%s.cutree.txt",out.prefix),sep = "\t",row.names = F,quote = F)
            ### make plot
            pdf(sprintf("%s.cutree.pdf",out.prefix),width=10,height=12)
            par(mar=c(5,4,4,2))
            layout(matrix(c(1,2),nrow = 2),heights = c(0.6,0.4))
            plot(obj.hclust.col,sub="",xlab="",hang=-1,cex=1.0*50/max(m,32))
            par(mar=c(2,4,0,2))
            colSet.cls <- auto.colSet(10)
            col.cls <- t(apply(dend.cutree,1,function(x){ colSet.cls[x] }))
            plotHclustColors(obj.hclust.col, colors=col.cls, cex.rowLabels = 1.1)
            dev.off()
        }
        pdf(sprintf("%s.pdf",out.prefix),width=pdf.width,height=pdf.height)
        par(mar=c(5,14,4,4))
        plot.new()
        title(main = mytitle,cex.main=4)
        #legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
        ### Integrating Grid Graphics Output with Base Graphics Output
        vps <- baseViewports()
        pushViewport(vps$inner, vps$figure, vps$plot)
        
        annDF <- data.frame(sampleType=sampleType)
        annColList <- list(sampleType=colSet)
        if(!is.null(ann.extra.df))
        {
            annDF <- cbind(annDF,ann.extra.df)
            annColList <- c(annColList,ann.extra.df.col)
        }
        if(!is.null(clonotype.col))
        {
            annDF$clonotype=clonotype.col$ctypeII[colnames(dat.plot)]
            annColList$clonotype=clonotype.col$colII
        }
        g.show.legend <- T
        if(do.cuttree){
            annDF <- cbind(annDF,dend.cutree)
            for(i in seq_len(ncol(dend.cutree))){
                tmp.list <- list(structure(colSet.cls,names=seq_along(colSet.cls)))
                names(tmp.list) <- colnames(dend.cutree)[i]
                annColList <- c(annColList,tmp.list)
                ### list(TReg=list(color <- bar="continuous",legend <- width=unit(2, "cm"),legend <- height=unit(4, "cm"))
                annotation_legend_param <- c(annotation_legend_param,list())
            }
            ann.bar.height <- max(0.5,ann.bar.height/ncol(annDF))
            g.show.legend <- F
            legend("topright",legend=names(colSet),fill=colSet,border=colSet,cex=1.5,inset=c(-0.03,0),xpd=T)
        }
        #print(str(annDF))
        #print(str(annColList))

        if(save.obj){
            save(dat.plot,z.title,z.lo,z.hi,z.step,do.scale,branch.col,branch.row,
                 pdf.width,pdf.height,mytitle,annDF,annColList,annotation_legend_param,ann.bar.height,
                 file=sprintf("%s.RData",out.prefix))
        }
        if(!is.null(sampleType) || !is.null(ann.extra.df)){
            ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = g.show.legend, annotation_legend_param = annotation_legend_param)
        }else{
            ha.col <- NULL
        }
        top_annotation_height <- unit(ann.bar.height * ncol(annDF), "cm")
        ###col = colorRamp2(c(0, 0.609, 1, 10), c("darkblue", "darkblue", "yellow", "red")),
        ####
                    ###col = colorRamp2(seq(-bk.range[2],bk.range[2],length=5), bluered(5),space="LAB"),
        bk.range <- quantile(abs(dat.plot),probs=c(0.01,1),na.rm=T)
        ht <- Heatmap(dat.plot, name=z.title,
                    ##col = colorRamp2(seq(-bk.range[2],bk.range[2],length=100), 
                    col = colorRamp2(seq(z.lo,z.hi,length=100), 
                                     colorRampPalette(rev(brewer.pal(n = 7, name = ifelse(do.scale,"RdBu","RdYlBu"))))(100), space="LAB"),
                                     ####colorRampPalette(rev(brewer.pal(n = 7, name =  "RdBu")))(100), space="LAB"),
                    column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"),
                    column_names_gp = gpar(fontsize = 12*55/max(m,32)),row_names_gp = gpar(fontsize = 10*55/max(n,32)),
                    show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
                    top_annotation_height = top_annotation_height,
                    cluster_columns = branch.col,
                    cluster_rows = branch.row,
                    heatmap_legend_param = list(grid_width = unit(0.8, "cm"), 
                                                grid_height = unit(0.8, "cm"), 
                                                at = seq(z.lo,z.hi,z.step),
                                                title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                label_gp = gpar(fontsize = 12), color_bar = "continuous"),
                    top_annotation = ha.col)
        ##print(str(draw))
        ##print((draw))
        ##print(class(ht))
        ComplexHeatmap::draw(ht, newpage= FALSE)
        for(i in seq_along(names(annColList))){
            decorate_annotation(names(annColList)[i], 
                                {grid.text(names(annColList)[i], unit(-4, "mm"),gp=gpar(fontsize=24),just = "right")})
        }
        dev.off()
    }
}

qcAndViz<-function(vsd,vstMat,designM,outDir,extra="",sfilter=NULL, gfilter=NULL,pair=FALSE)
{
    if(!is.null(sfilter))
    {
	    vstMat <- vstMat[,sfilter]
	    designM <- designM[sfilter,]
    }
    if(!is.null(gfilter))
    {
	    vstMat <- vstMat[rownames(vstMat) %in% gfilter[!is.na(gfilter)],]
    }
	vstMat <- vstMat[, rownames(designM)]
#    ## visualization
#    pdf(paste(outDir,"/DESeq2.plot.QA.pdf",sep=""))
#    ## MA plot
#    plotMA(res, main="DESeq2", ylim=c(-2,2))
#    ## dispersion plot
#    plotDispEsts(x)
#    dev.off()
    ## PCA
    suppressPackageStartupMessages(require("DESeq2"))
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    cat(sprintf("%s\t---- PCA plot\n", Sys.time()))
    if(pair)
    { 
        plotPCA.dat <- plotPCA(vsd, intgroup=c("sampleType", "patient"),returnData=T)
	    percentVar <- round(100 * attr(plotPCA.dat, "percentVar"))
            ggplot(data = plotPCA.dat, aes_string(x = "PC1", y = "PC2", color = "group")) + geom_point(size = 3) + 
		    xlab(paste0("PC1: ", round(percentVar[1]), "% variance")) + 
		    ylab(paste0("PC2: ", round(percentVar[2]), "% variance"))
	    ggsave(filename=paste(outDir,"/DESeq2.PCA",extra,".pdf",sep=""),width=8,height=6)
	    write.table(plotPCA.dat,paste(outDir,"/DESeq2.PCA",extra,".txt",sep=""),col.names = T,row.names = F,quote = F)

    }else
    {
        plotPCA.dat <- plotPCA(vsd, intgroup=c("sampleType"),returnData=T) 
	    percentVar <- round(100 * attr(plotPCA.dat, "percentVar"))
            ggplot(data = plotPCA.dat, aes_string(x = "PC1", y = "PC2", color = "group")) + geom_point(size = 3) + 
		    xlab(paste0("PC1: ", round(percentVar[1]), "% variance")) + 
		    ylab(paste0("PC2: ", round(percentVar[2]), "% variance"))
	    ggsave(filename=paste(outDir,"/DESeq2.PCA",extra,".pdf",sep=""),width=8,height=6)
	    write.table(plotPCA.dat,paste(outDir,"/DESeq2.PCA",extra,".txt",sep=""),col.names = T,row.names = F,quote = F)

    }
    ## heatmap by most variable genes
    cat(sprintf("%s\t---- clustering variable genes\n", Sys.time()))

    if(length(levels(designM$sampleType))<=9)
    {
	    suppressPackageStartupMessages(require("RColorBrewer"))
	    colSet <- brewer.pal(9,"Set1")
	    patientcolors <- colSet[as.numeric(designM$sampleType)]
    }else
    {
        colSet <- rainbow(length(levels(designM$sampleType)))
	    patientcolors <- colSet[as.numeric(designM$sampleType)]
    }
    #print(patientcolors)
    designM$color <- patientcolors
    legendDat <- unique(designM[,c("sampleType","color")])
    
    rowVar <- apply(vstMat,1,var)
    for(n in c(50,100,150,200,250,300,350,400,450,500))
    {
        ###hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(n, length(rowVar)))]
        dat.plot <-vstMat[select,]
        
        dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=entrezToXXX(rownames(dat.plot)))
        dat.plot.export <- cbind(dat.plot.export,dat.plot)
	    write.table(dat.plot.export,paste(outDir,"/DESeq2.QC.Cluster.Var.n",n,extra,".txt",sep=""),sep="\t",col.names = T,row.names = F,quote = F)
        
        nn <- min(n, length(rowVar))
        m <- ncol(dat.plot)
        rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
        pdf(paste(outDir,"/DESeq2.QC.Cluster.Var.n",n,extra,".pdf",sep=""),width=12,height=12)
        heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20), main="Most variable genes",cexRow=min(1.8,55/nn),cexCol=min(1.8,55/m))
        legend("topright",legend=legendDat$sampleType,fill=legendDat$color,border=legendDat$color,cex=1.5)
        dev.off()
    }
    cat(sprintf("%s\t---- qcAndViz finished\n", Sys.time()))

}

qcAndVizMat<-function(vstMat,designM,outDir,colSet=NULL,intgroup=c("sampleType"),ntop=500,extra="",sfilter=NULL, gfilter=NULL,complexHeatmap.use=NULL,clonotype.col=NULL,runNMF=FALSE,...)
{
    dir.create(outDir,showWarnings = F,recursive = T)
    if(!is.null(sfilter))
    {
	    vstMat <- vstMat[, colnames(vstMat) %in% sfilter, drop = FALSE]
	    designM <- designM[rownames(designM) %in% sfilter,, drop = FALSE]
    }
    if(!is.null(gfilter))
    {
	    vstMat <- vstMat[rownames(vstMat) %in% gfilter[!is.na(gfilter)],, drop = FALSE]
    }
    ###print(dim(designM))
    ###print(dim(vstMat))
	vstMat <- vstMat[, rownames(designM), drop = FALSE]
    ###print(str(vstMat))
    ###return(NULL)

    if(is.null(colSet))
    {
        nLevel <- length(levels(designM$sampleType))
        if(nLevel<=9)
        {
            suppressPackageStartupMessages(require("RColorBrewer"))
            colSet <- brewer.pal(nLevel,"Set1")
            names(colSet) <- levels(designM$sampleType)
        }else
        {
            colSet <- auto.colSet(n=nLevel)
            names(colSet) <- levels(designM$sampleType)
        }
    }else
    {
      colSet <- colSet[names(colSet) %in% levels(designM$sampleType)]
    }
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("factoextra"))
    rowVar <- apply(vstMat,1,var)
    f.var <- rowVar > 0
    vstMat <- vstMat[f.var,]
    rowVar <- rowVar[f.var]
    ###print(dim(vstMat))
    #rowVar <- apply(vstMat,1,function(x){ var(x)/(mean(x)^2) } )
    if(is.null(ntop))
    {
        ntop=length(rowVar)
    }
    select <- order(rowVar, decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
    
    ## PCA
    runPCAAnalysis(vstMat[select, ],sprintf("%s/PCA%s",outDir,extra),
                   designM$sampleType,colSet,ntop=NULL,main=sub("^.","",extra))
    ## t-SNE
    runTSNEAnalysis(vstMat[select, ],sprintf("%s/tsne%s",outDir,extra),
                    names(colSet),colSet[as.character(designM$sampleType)],colSet)
    ### hclustering
    for(n in unique(c(50,100,150,200,250,500,ntop)))
    #for(n in c(50,100,150,200,250,300,350,400,450,500,1000,2000,ntop))
    {
        if(!is.null(clonotype.col))
        {
            runHierarchicalClusteringAnalysis(vstMat,
                                              sprintf("%s/QC.Cluster.Var.n%s%s.share",outDir,n,extra),
                                              designM$sampleType,colSet,
                                              clonotype.col=clonotype.col[["share"]],
                                              ntop=n,
                                              complexHeatmap.use=complexHeatmap.use,
                                              verbose=TRUE,main="Variable Genes")
            runHierarchicalClusteringAnalysis(vstMat,
                                              sprintf("%s/QC.Cluster.Var.n%s%s.strict",outDir,n,extra),
                                              designM$sampleType,colSet,
                                              clonotype.col=clonotype.col[["strict"]],
                                              ntop=n,
                                              complexHeatmap.use=complexHeatmap.use,
                                              verbose=TRUE,main="Variable Genes")
        }else{
            runHierarchicalClusteringAnalysis(vstMat,
                                              sprintf("%s/QC.Cluster.Var.n%s%s",outDir,n,extra),
                                              designM$sampleType,colSet,
                                              clonotype.col=clonotype.col,
                                              ntop=n,
                                              complexHeatmap.use=complexHeatmap.use,
                                              verbose=TRUE,main="Variable Genes")
        }
        ## sample distance
#        if(!file.exists(paste(outDir,"/DESeq2.QC.SampleDist",extra,".pdf")))
#        {
#            distsRL <- dist(t(vstMat))
#            mat <- as.matrix(distsRL)
#            pdf(paste(outDir,"/DESeq2.QC.SampleDist",extra,".pdf",sep=""),width=12,height=12)
#            heatmap.2(mat, RowSideColors=patientcolors, ColSideColors=patientcolors, density.info="none", keysize=1.2, trace="none", margin=c(15, 15), main="Sample Distance",cexRow=1.8,cexCol=1.8)
#            dev.off()
#        }
    }
    # NMF
    if(runNMF)
    {
        runNMFAnalysis(vstMat[select, ],sprintf("%s/NMF%s",outDir,extra),
                       designM[,"sampleType",drop=F],list(sampleType=colSet))
    }
    loginfo("qcAndVizMat() finished")
}


#read.ExpData <- function(exp.file,design.file) {
read.ExpData <- function(exp.file,design=NULL) {
  if( grepl(".RData$",tab.file,perl=T) )
  {
      suppressPackageStartupMessages(library("R.utils"))
      lenv <- loadToEnv(exp.file)
      #vsd.blind <- lenv[["vsd.blind"]]
      in.table <- lenv[["vstMat.blind"]]

  }else{
      in.table <- read.table(exp.file,sep="\t",stringsAsFactors = F,header = T,row.names = 1,check.names = F)
      in.table <- as.matrix(in.table[,-1])
  }
  #design<-read.table(design.file,header=T,row.names="sample",check.names=F,colClasses=c("factor","character","factor","factor"))
  if(is.null(design))
  {
      in.table
  }else
  {
      in.table[,rownames(design)]
  }
}

read.SampleTypeColor <- function(in.file) {
    if(is.null(in.file) || !file.exists(in.file)) { return(NULL); }
    in.table <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F,sep="\t",comment.char="!")
    ret <- in.table$color
    names(ret) <- in.table$sampleType
    return(ret)
}

patientColorListFromMyDesign <- function(d) {
    ctype <- d$libType
    names(ctype) <- rownames(d)
    c.palette.fct <- as.factor(d$libType)
    c.palette <- seq_along(levels(c.palette.fct))+1
    names(c.palette) <- levels(c.palette.fct)
    list(ctype=ctype,col=c.palette)
}

read.clonotype <- function(in.file,ctype.col) {
    suppressPackageStartupMessages(require("RColorBrewer"))
    if(is.null(in.file) || !file.exists(in.file)) { return(NULL) }
    in.table <- read.table(in.file,row.names = "Cell_Name",header = T,check.names = F,stringsAsFactors = F,sep="\t",comment.char="!",quote = "")
    ctype <- in.table[,ctype.col]
    ctype <- strsplit(x = ctype,split = ":",perl = T)
    ###ctype <- sapply(ctype,function(x){ if(as.integer(x[2])>1) { x[1] } else { "NoClonal" }  })
    ctypeI <- sapply(ctype,function(x){ if(as.integer(x[2])>1) { "Clonal" } else { "NoClonal" }  })
    ctypeII <- sapply(ctype,function(x){ n=as.integer(x[2]); if(n>4) n=4; sprintf("N%d",n) })
    names(ctypeI) <- rownames(in.table)
    names(ctypeII) <- rownames(in.table)
    #ctypeI <- as.factor(ctypeI)
    #ctypeII <- as.factor(ctypeII)
    c.type.level.I <- unique(ctypeI)
    c.type.level.II <- unique(ctypeII)
    ##c.palette <- rainbow(length(c.type.level))
    #names(c.palette) <- c.type.level
    ##c.palette[names(c.palette)=="NoClonal"]="gray60"
    ##c.palette[names(c.palette)=="Clonal"]="#7FC97F"
    c.palette.I <- structure(auto.colSet(n = length(c.type.level.I) ,name = "Paired"),names=c.type.level.I)
    ###c.palette.II <- structure(auto.colSet(n = length(c.type.level.II) ,name = "Paired"),names=c.type.level.II)
    c.palette.II <- structure(brewer.pal(11,"BrBG")[c(7,8,9,10)],names=c("N1","N2","N3","N4"))
    c.palette.I <- c(c.palette.I, structure(c("gray"), names = c("NA")))
    c.palette.II <- c(c.palette.II, structure(c("gray"), names = c("NA")))
    #c.color <- c.palette.I[ctypeI]
    #names(c.color) <- rownames(in.table)
    list(ctype=ctypeI,col=c.palette.I,ctypeII=ctypeII,colII=c.palette.II)
}

runTopGOAnalysis <- function( goiIDs, universeIDs , mapping.db="org.Hs.eg.db" )
{
    require("topGO")
    sapply( c( "MF", "BP", "CC" ), function( ont ) { 
           alg <- factor( as.integer( universeIDs %in% goiIDs ) )
           print(str(alg))
           names(alg) <- goiIDs
           tgd <- new( "topGOdata", ontology=ont, allGenes = alg, nodeSize=5, annot=annFUN.org, mapping=mapping.db)
           resultTopGO <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
           GenTable( tgd, resultTopGO, topNodes=15 ) 
       }, simplify=FALSE )
}

require("methods")
require("DESeq2")
## class definition
SCDenoise <- setClass("SCDenoise", slots = c(raw.endo = "matrix", normalized.endo = "matrix", final.endo = "matrix", size.factor.endo = "vector",
                                             raw.ERCC = "matrix", normalized.ERCC = "matrix", final.ERCC = "matrix", size.factor.ERCC = "vector",
                                             withERCC = "logical"))
setValidity("SCDenoise",
            function(object) {
              msg <- NULL
              if ( ! is.matrix(object@raw.endo) ){
                msg <- c(msg, "input data must be data.frame")
              }else if ( nrow(object@raw.endo) < 2 ){
                msg <- c(msg, "input data must have more than one row")
              }else if ( ncol(object@raw.endo) < 2 ){
                msg <- c(msg, "input data must have more than one column")
              }
              if (is.null(msg)) TRUE
              else msg
            }
            )
setMethod("initialize",
          signature = "SCDenoise",
          definition = function(.Object, expdata, ignore.ERCC=FALSE){
            ###.Object@normalized.endo <- expdata
            ###.Object@final.endo <- expdata
            .Object@raw.ERCC <- expdata[grepl(pattern = "^ERCC-",x = rownames(expdata),perl = T),]
            .Object@raw.endo <- expdata[!grepl(pattern = "^ERCC-",x = rownames(expdata),perl = T),]
            if(nrow(.Object@raw.ERCC)==0 || ignore.ERCC==TRUE) { 
                .Object@withERCC <- FALSE 
            }else { 
                .Object@withERCC <- TRUE 
            }
            validObject(.Object)
            return(.Object)
          }
          )

setGeneric("SCDenoise.normalize", function(object,useERCCSizeFactor=FALSE) standardGeneric("SCDenoise.normalize"))
setMethod("SCDenoise.normalize",
          signature = "SCDenoise",
          definition = function(object,useERCCSizeFactor=FALSE) {
              endo.sf <- estimateSizeFactorsForMatrix(object@raw.endo)
              object@size.factor.endo <- endo.sf
              if(object@withERCC==FALSE) { 
                  useERCCSizeFactor <- FALSE 
                  ######object@size.factor.ERCC <- NULL
                  object@normalized.ERCC <- object@raw.ERCC
              }else
              {
                  try({
                      ERCC.sf <- estimateSizeFactorsForMatrix(object@raw.ERCC)
                      object@size.factor.ERCC <- ERCC.sf
                      object@normalized.ERCC <- t( t(object@raw.ERCC) / ERCC.sf )
                  })
              }
              if(useERCCSizeFactor) {
                  object@normalized.endo <- t( t(object@raw.endo) / ERCC.sf )
              }else{
                  object@normalized.endo <- t( t(object@raw.endo) / endo.sf )
              }
              return(object)
          }
          )

setGeneric("SCDenoise.fitTechnicalNoise", function(object, use_ERCC = TRUE,fit_type = "counts",plot=TRUE, fit_opts=NULL) standardGeneric("SCDenoise.fitTechnicalNoise"))
setMethod("SCDenoise.fitTechnicalNoise",
          signature = "SCDenoise",
          definition = function(object, use_ERCC = TRUE,fit_type = "counts",plot=TRUE, fit_opts=NULL) {
                nCountsEndo <- object@normalized.endo
                nCountsERCC <- object@normalized.ERCC
                #### check for parameters
                if(object@withERCC==FALSE && use_ERCC==TRUE) { 
                    cat("No ERCC available, use endo-genes to do fitting (set use_ERCC to FALSE)\n")
                    use_ERCC <- FALSE 
                }
                ##if( use_ERCC==FALSE &&  (fit_type %in% c('counts','logvar'))){
                ##    warning("Without ERCCs 'fit_type' 'log' is recommedned")
                ##}
                if((fit_type %in% c('counts', 'log','logvar'))==F){stop("'fit_type' needs to be 'counts', 'log' or 'logvar'")}
                ##if(fit_type=="counts" & use_ERCC==FALSE){
                ##   print("Without ERCCs fit needs to be perfromed in log-space")
                ##   use_ERCC = FALSE  
                ##}
 
                if(use_ERCC==TRUE){
                    if(fit_type=="counts"){
                        meansEndo <- rowMeans( nCountsEndo )
                        varsEndo <- rowVars( nCountsEndo )
                        cv2Endo <- varsEndo / meansEndo^2

                        meansERCC <- rowMeans( nCountsERCC )
                        varsERCC <- rowVars( nCountsERCC )
                        cv2ERCC <- varsERCC / meansERCC^2
      
                        #Do fitting of technical noise
                        if(!is.null(fit_opts)){
                            if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
                            if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
                        }else{
                            mincv2 = 0.3
                            quan = 0.8
                        }
      
                        #normalised counts (with size factor)
                        minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > mincv2 ) ], quan ) )
                        useForFitA <- meansERCC >= minMeanForFitA
                        fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ), cv2ERCC[useForFitA] )
     
                        ###xi <- mean( 1 / object@size.factor.ERCC )
                        ###a0 <- unname( fitA$coefficients["a0"] )
                        ###a1 <- unname( fitA$coefficients["a1tilde"] - xi ) 
                        #4. Transform to log-space and propagate error
                        eps=1
                        LogNcountsEndo=log10(nCountsEndo+eps)
                        dLogNcountsEndo=1/((meansEndo+eps)*log(10))
                        var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
                        LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
      
                        if(plot==TRUE){
                            #plot fit
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( meansERCC, cv2ERCC, log="xy", col=1+2*useForFitA, pch=19, xlab = 'Means', ylab = expression("CV" ^ 2))
                            xg <- 10^seq( -3, 5, length.out=100 )
                            lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='blue' )
                            segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
                            meansERCC[useForFitA], fitA$fitted.values, col="gray" )
                            legend('bottomleft',c('Genes used for fit', 'Fit technical noise'),pch=c(19, NA),lty =c(NA,1),col=c('green','blue'),cex=0.8)
                            title(expression("CV" ^ 2 * " ~ Mean (using ERCC)"))
                            
                            #plot fot with all genes
                            plot( meansEndo, cv2Endo, log="xy", col=1, xlab = 'Means', ylab = expression("CV" ^ 2), pch=20, cex=0.3)
                            points(meansERCC, cv2ERCC, col='#0000FFF0', pch=15, cex=0.6)
                            xg <- 10^seq( -3, 5, length.out=100 )
                            lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='#0000FFF0',lwd=2 )

                            # Plot quantile lines around the fit
                            df <- ncol(nCountsEndo) - 1
                            lines( xg, ( coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg ) * qchisq( .975, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )
                            lines( xg, ( coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg ) * qchisq( .025, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )

                            legend('bottomleft',c('Endogenous genes','ERCCs', 'Fit technical noise'),pch=c(20,15, NA),lty =c(NA,NA,1),col=c('black','blue','blue'),cex=0.8)        
                            title(expression("CV" ^ 2 * " ~ Mean (using ERCC)"))
                            ###par(mfrow=c(1,1))
                        }
                        res = list()
                        res$fit = fitA
                        res$techNoiseLog = LogVar_techEndo
                    }else{#with ERCCs in log space
                        if(fit_type=="log"){
                            LCountsEndo <- log10(nCountsEndo+1)
                            LmeansEndo <- rowMeans( LCountsEndo )
                            LvarsEndo <- rowVars( LCountsEndo )
                            Lcv2Endo <- LvarsEndo / LmeansEndo^2
                            
                            LCountsERCC = log10(nCountsERCC+1)
                            LmeansERCC <- rowMeans( LCountsERCC )
                            LvarsERCC <- rowVars( LCountsERCC )
                            Lcv2ERCC <- LvarsERCC / LmeansERCC^2
                            
                            if(!is.null(fit_opts)){
                                if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=2}
                            }else{
                                minmean = .5
                            }
                            LogNcountsList=list()
                            useForFitL=LmeansERCC>minmean
                            LogNcountsList$mean=LmeansERCC[useForFitL]
                            LogNcountsList$cv2=Lcv2ERCC[useForFitL]
                            fit_loglin=nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=20,k=1))
                            LogVar_techEndo_logfit <- coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansEndo)*LmeansEndo^2
                            
                            if(plot==TRUE){
                                plot( LmeansEndo, Lcv2Endo, log="y", col=1,ylim=c(1e-3,1e2),xlab='meansLogEndo',ylab='cv2LogEndo')
                                xg <- seq( 0, 5.5, length.out=100 )
                                lines( xg, coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),lwd=2,col='green' )
                                points(LmeansERCC, Lcv2ERCC,col='blue',pch=15,cex=1.1)
                                legend('topright',c('Endo','ERCC'),pch=c(1,1,15),col=c('black','blue'))
                            }
                            res = list()
                            res$fit = fit_loglin
                            res$techNoiseLog = LogVar_techEndo_logfit
                        }else{#with ERCCs fit variance in log space with loess
                            LCountsEndo <- log10(nCountsEndo+1)
                            LmeansEndo <- rowMeans( LCountsEndo )
                            LvarsEndo <- rowVars( LCountsEndo )
                            Lcv2Endo <- LvarsEndo / LmeansEndo^2
                            
                            LCountsERCC = log10(nCountsERCC+1)
                            LmeansERCC <- rowMeans( LCountsERCC )
                            LvarsERCC <- rowVars( LCountsERCC )
                            
                            if("span" %in% names(fit_opts)){span = fit_opts$span}else{span=0.8}
                            if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=0.5}
                            
                            useForFitA <- LmeansERCC >= minmean
                            fit_var2 = loess(LvarsERCC[useForFitA] ~ LmeansERCC[useForFitA], span=span, control=loess.control(surface="direct"))
                            xg <- seq( 0, 5.5, length.out=100 )
                            Var_techEndo_logfit_loess <-  predict(fit_var2, xg)
                            
                            minVar_ERCC = min(LvarsERCC[LmeansERCC>3])
        
                            if(any(xg>3 & (Var_techEndo_logfit_loess<0.8*minVar_ERCC))){
                                idx_1 = which(xg>3 & (Var_techEndo_logfit_loess<0.8*minVar_ERCC))[1]
                                idx_end = length(Var_techEndo_logfit_loess)
                                Var_techEndo_logfit_loess[idx_1:idx_end] = 0.8*minVar_ERCC        
                            }
                            
                            if(plot==TRUE){
                                plot( LmeansEndo, LvarsEndo, col=1,ylim=c(1e-3,150.5),log="y",xlab='meansLogEndo',ylab='VarLogEndo')
                                points(LmeansERCC, LvarsERCC,col='blue',pch=15,cex=1.1)
                                lines(xg, Var_techEndo_logfit_loess,lwd=3,col='blue',lty=1)  
                                legend('topright',c('Endo. genes','ERCC', 'Tech. noise fit'),pch=c(1,15,NA), lty = c(NA,NA,1),col=c('black','blue', 'blue'))
                            }
                            
                            #use model for endogenous genes
                            xg=LmeansEndo
                            Var_techEndo_logfit_loess <-  predict(fit_var2, xg)      
                            
                            if(any(xg>3 & Var_techEndo_logfit_loess<0.8*minVar_ERCC)){
                            idx_1 = which(xg>3 & Var_techEndo_logfit_loess<0.8*minVar_ERCC)[1]
                            idx_end = length(Var_techEndo_logfit_loess)
                            Var_techEndo_logfit_loess[idx_1:idx_end] = 0.8*minVar_ERCC       
                            }          
                            
                            res = list()
                            res$fit = fit_var2
                            res$techNoiseLog = Var_techEndo_logfit_loess
                        }
                    }
                }else{#no ERCCs available
                    if(fit_type=="log"){
                        LCountsEndo <- log10(nCountsEndo+1)
                        LmeansEndo <- rowMeans( LCountsEndo )
                        LvarsEndo <- rowVars( LCountsEndo )
                        Lcv2Endo <- LvarsEndo / LmeansEndo^2
                        
                        if(!is.null(fit_opts)){
                            if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{fit_opts$minmean=0.3}
                            if("offset" %in% names(fit_opts)){offset = fit_opts$offset}else{fit_opts$offset=1}
                        }else{
                            fit_opts$minmean = 0.3
                            fit_opts$offset=1
                        }
                        
                        LogNcountsList = list()
                        useForFitL = LmeansEndo>fit_opts$minmean
                        LogNcountsList$mean = LmeansEndo[useForFitL]
                        LogNcountsList$cv2 = Lcv2Endo[useForFitL]
                        fit_loglin = nls(cv2 ~ a* 10^(-k*mean), LogNcountsList,start=c(a=10,k=2))
                        fit_loglin$opts = fit_opts      
                        LogVar_techEndo_logfit <- fit_opts$offset* coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*LmeansEndo)*LmeansEndo^2
                        
                        if(plot==TRUE){
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( LmeansEndo, Lcv2Endo, log="y", col=1,ylim=c(1e-3,1e2),
                                 xlab='means',ylab=expression("CV" ^ 2),
                                 main=expression("CV" ^ 2 * " ~ Mean (using Endo in log space)"),pch=20,cex=0.3)
                            xg <- seq( 0, 5.5, length.out=100 )
                            lines( xg, fit_opts$offset*coefficients(fit_loglin)["a"] *10^(-coefficients(fit_loglin)["k"]*xg ),
                                  lwd=2,col='green' )
                        }
                        res = list()
                        res$fit = fit_loglin
                        res$techNoiseLog = LogVar_techEndo_logfit
                    }
                    if(fit_type=='counts'){
                        meansEndo <- rowMeans( nCountsEndo )
                        varsEndo <- rowVars( nCountsEndo )
                        cv2Endo <- varsEndo / meansEndo^2
                        
                        #Do fitting of technical noise
                        if(!is.null(fit_opts)){
                            if("mincv2" %in% names(fit_opts)){mincv2 = fit_opts$mincv2}else{mincv2=.3}
                            if("quan" %in% names(fit_opts)){quan = fit_opts$quan}else{quan=0.8}
                        }else{
                            mincv2 = 0.3
                            quan=0.8
                        }
                        
                        #normalised counts (with size factor)
                        minMeanForFitA <- unname( quantile( meansEndo[ which( cv2Endo > mincv2 ) ], quan ) )
                        useForFitA <- meansEndo >= minMeanForFitA
                        fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansEndo[useForFitA] ),
                        cv2Endo[useForFitA] )
                        
                        #4. Transform to log-space and propagate error
                        eps=1
                        LogNcountsEndo=log10(nCountsEndo+eps)
                        dLogNcountsEndo=1/((meansEndo+eps)*log(10))
                        var_techEndo=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansEndo)*meansEndo^2
                        LogVar_techEndo=(dLogNcountsEndo*sqrt(var_techEndo))^2 #error propagation 
                        
                        if(plot==TRUE){
                            #plot fit        
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( meansEndo, cv2Endo, log="xy", col=1+2*useForFitA, pch=20, cex=0.3, xlab = 'Means', ylab = expression("CV" ^ 2))
                            xg <- 10^seq( -3, 5, length.out=100 )
                            lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg, col='#0000FFF0' )
                            legend('bottomleft',c('Genes used for fit', 'Fit baseline variation'),
                                   pch=c(20, NA),lty =c(NA,1),col=c('green','blue'),cex=1.2)
                            title(expression("CV" ^ 2 * " ~ Mean (using endogeneous genes)"))
                        }
                        res = list()
                        res$fit = fitA
                        res$techNoiseLog = LogVar_techEndo
                    }
                    if(fit_type=='logvar'){
                        LCountsEndo <- log10(nCountsEndo+1)
                        LmeansEndo <- rowMeans( LCountsEndo )
                        LvarsEndo <- rowVars( LCountsEndo )
                        Lcv2Endo <- LvarsEndo / LmeansEndo^2
                        
                        if("span" %in% names(fit_opts)){span = fit_opts$span}else{span=0.8}
                        if("minmean" %in% names(fit_opts)){minmean = fit_opts$minmean}else{minmean=0.5}
                        
                        useForFitA <- LmeansEndo >= minmean
                        fit_var2 = loess(LvarsEndo[useForFitA] ~ LmeansEndo[useForFitA], span=span, control=loess.control(surface="direct"))
                        xg <- seq( 0, 5.5, length.out=100 )
                        Var_techEndo_logfit_loess <-  predict(fit_var2, xg)
                        
                        minVar_ERCC = min(LvarsEndo[LmeansEndo>3])
                        
                        if(any(xg>2.5 & (Var_techEndo_logfit_loess<0.6*minVar_ERCC))){
                            idx_1 = which(xg>2.5 & (Var_techEndo_logfit_loess<0.6*minVar_ERCC))[1]
                            idx_end = length(Var_techEndo_logfit_loess)
                            Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_ERCC        
                        }
                        
                        if(plot==TRUE){
                            par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                            plot( LmeansEndo, LvarsEndo, col=1,ylim=c(1e-3,150.5),log="y",xlab='means',ylab='Var',main=expression("Var~ Mean (using Endo in log space)"),pch=20,cex=0.3)
                            lines(xg, Var_techEndo_logfit_loess,lwd=3,col='blue',lty=1)  
                            legend('topright',c('Endo. genes', 'Tech. noise fit'),pch=c(1,NA), lty = c(NA,1),col=c('black', 'blue'))
                        }
                        
                        #use model for endogenous genes
                        xg=LmeansEndo
                        Var_techEndo_logfit_loess <-  predict(fit_var2, xg)      
                        
                        if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_ERCC)){
                            idx_1 = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_ERCC)[1]
                            idx_end = length(Var_techEndo_logfit_loess)
                            Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_ERCC       
                        }  
                        
                        res = list()
                        res$fit = fit_var2
                        res$techNoiseLog = Var_techEndo_logfit_loess
                    }
                }
                res$fit_opts = fit_opts
                res  
            })

setGeneric("SCDenoise.getVariableGenes", function(object, fit, method = "fit", threshold = 0.1, fit_type=NULL, plot=T, fitB=NULL) standardGeneric("SCDenoise.getVariableGenes"))
setMethod("SCDenoise.getVariableGenes",
            signature = "SCDenoise",
            definition = function(object, fit, method = "fit", threshold = 0.1, fit_type=NULL, plot=T, fitB=NULL) {
                nCountsEndo <- object@normalized.endo
                nCountsERCC <- object@normalized.ERCC
                sfEndo <- object@size.factor.endo
                sfERCC <- object@size.factor.ERCC
                res.list <- list()
                if(!(method %in% c("fdr","fit"))){ stop("'method' needs to be either 'fdr' or 'fit'") }
                if(is.null(fit_type)){
                  print("No 'fit_type' specified. Trying to guess its from parameter names")
                  if("a0" %in% names(coefficients(fit)) & "a1tilde" %in% names(coefficients(fit))){fit_type="counts"}else{
                    if("a" %in% names(coefficients(fit)) & "k" %in% names(coefficients(fit))){fit_type="log"}else{
                      if(is.call(fit$call)){fit_type="logvar"}
                    }
                  }
                  print(paste("Assuming 'fit_type' is ","'",fit_type,"'",sep=""))
                }
                if(is.null(fit_type)){stop("Couldn't guess fit_type. Please specify it or run the fitTechnicalNoise function to obtain the fit")}
                if(!(fit_type %in% c("counts","log", "logvar")) & !is.null(fit_type)){ stop("'fit_type' needs to be either 'fdr' or 'fit'") }
                if(method=='fdr' & fit_type!="counts"){stop("method='fdr', can only be used with fit_type 'counts'")}
                if(method=='fdr' & object@withERCC==FALSE) { method <- "fit" }
                if(method=='fdr' & (is.null(sfERCC) | is.null(sfEndo))){stop("Please specify sfERCC and sfEndo when using method='fdr'")}
                if(method=='fdr'){
                  meansEndo <- rowMeans( nCountsEndo )
                  varsEndo <- rowVars( nCountsEndo )
                  cv2Endo <- varsEndo / meansEndo^2
                   
                  meansERCC <- rowMeans( nCountsERCC )
                  varsERCC <- rowVars( nCountsERCC )
                  cv2ERCC <- varsERCC / meansERCC^2
              
                  minBiolDisp <- .5^2
                    
                  xi <- mean( 1 / sfERCC )
                  ###a0 <- unname( fit$coefficients["a0"] )
                  ###a1 <- unname( fit$coefficients["a1tilde"] - xi )
                  m <- ncol(nCountsEndo)
                  psia1thetaA <- mean( 1 / sfEndo ) + ( coefficients(fit)["a1tilde"] - xi ) * mean( sfERCC / sfEndo )
                  cv2thA <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
                  testDenomA <- ( meansEndo * psia1thetaA + meansEndo^2 * cv2thA ) / ( 1 + cv2thA/m )
                  pA <- 1 - pchisq( varsEndo * (m-1) / testDenomA, m-1 )
                  padjA <- p.adjust( pA, "BH" )
                  #print( table( padjA < .1 ))
                  is_het =  padjA < threshold
                  is_het[is.na(is_het)] = FALSE
                  res.list[["padjA"]] <- padjA
                  res.list[["residual"]] <- cv2Endo - (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo)
                  res.list[["x"]] <- meansEndo
                  res.list[["y"]] <- cv2Endo
                  
                  if(plot==TRUE){
                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,250), xlab='Means', ylab = expression("CV" ^ 2), pch=20, cex=0.3)
                    xg <- 10^seq( -3, 5, length.out=100 )
                    lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='#0000FFF0' )
                    try({
                        points( meansERCC, cv2ERCC, col="blue", pch=15, cex=0.6 )
                        # Plot quantile lines around the fit
                        df <- ncol(nCountsEndo) - 1
                        lines( xg, ( coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg ) * qchisq( .975, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )
                        lines( xg, ( coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg ) * qchisq( .025, df ) / df, col="#0000FFF0", lwd=2, lty="dashed" )
                        lines( xg, psia1thetaA/xg + coefficients(fit)["a0"] + minBiolDisp, lty="dashed", col="#FF0000F0", lwd=2 )
                    })
                    if(!is.null(fitB))
                    {
                        print(rbind(byERCC=coefficients(fit),byEndo=coefficients(fitB)))
                        lines( xg, coefficients(fitB)["a0"] + coefficients(fitB)["a1tilde"]/xg, col='#00FF00F0', lwd=2 )
                        legend('bottomleft',c('Endo. genes','Var. genes','ERCCs',"Fit","Fit(Endo)"),
                               pch=c(20,20,15,NA,NA),lty = c(NA,NA,NA,1,1),
                               col=c('black','red','blue', '#0000FFF0', '#00FF00F0'),cex=1.2)
                    }else
                    {
                        legend('bottomleft',c('Endo. genes','Var. genes','ERCCs',"Fit"),
                               pch=c(20,20,15,NA),lty = c(NA,NA,NA,1),
                               col=c('black','red','blue', '#0000FFF0'),cex=1.2)
                    }
                    title(expression("Variable genes by Mean ~ CV" ^ 2 * " (using ERCCs)"))
                  }
                }
                if(method=='fit' & fit_type=='log')
                {
                  LCountsEndo <- log10(nCountsEndo+1)
                  LmeansEndo <- rowMeans( LCountsEndo )
                  Lcv2Endo = rowVars(LCountsEndo)/LmeansEndo^2
                  is_het = (fit$opts$offset * coefficients(fit)["a"] *10^(-coefficients(fit)["k"]*LmeansEndo) < Lcv2Endo) &  LmeansEndo>fit$opts$minmean 
                  res.list[["residual"]] <- Lcv2Endo - (fit$opts$offset * coefficients(fit)["a"] *10^(-coefficients(fit)["k"]*LmeansEndo))
                  res.list[["x"]] <- LmeansEndo
                  res.list[["y"]] <- Lcv2Endo
                  
                  if(plot==TRUE){
                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( LmeansEndo, Lcv2Endo, log="y", col=1+is_het,pch=20,cex=0.3,
                         xlab='means',ylab=expression("CV" ^2),main=expression("CV"^2*" ~ Means (log transformed)"))
                    xg <- seq( 0, 5.5, length.out=100 )
                    lines( xg, fit$opts$offset * coefficients(fit)[1] *10^(-coefficients(fit)[2]*xg ),lwd=2,col='green' )
                    legend('bottomright',c('Endo. genes','Var. genes',"Fit"),
                           pch=c(20,20,NA),lty = c(NA,NA,1),col=c('black','red', 'blue'),cex=1.2)
                  }
                }
                if(method=='fit' & fit_type=='counts')
                {
                  meansEndo <- rowMeans( nCountsEndo )
                  varsEndo <- rowVars( nCountsEndo )
                  cv2Endo <- varsEndo/meansEndo^2
                  is_het = (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo) < cv2Endo #&  meansEndo>2
                  res.list[["residual"]] <- cv2Endo - (coefficients(fit)[[1]] + coefficients(fit)[[2]]/meansEndo)
                  res.list[["x"]] <- meansEndo
                  res.list[["y"]] <- cv2Endo
                  
                  if(plot==TRUE){
                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( meansEndo, cv2Endo, log="xy", col=1+is_het,ylim=c(0.1,95), xlab='Means', ylab=expression("CV" ^ 2), pch=20, cex=0.3)
                    xg <- 10^seq( -3, 5, length.out=100 )
                    lines( xg, coefficients(fit)[1] + coefficients(fit)[2]/xg,lwd=2,col='#00FF00F0' )
                    legend('bottomright',c('Endo. genes','Var. genes',"Fit"),
                           pch=c(20,20,NA),lty = c(NA,NA,1),col=c('black','red', '#00FF00F0'),cex=1.2)   
                    title(expression("Variable genes by Mean ~ CV" ^ 2 * " (using endogenous genes)"))
                  }
                }
                if(method=='fit' & fit_type=='logvar')
                {
                  LCountsEndo <- log10(nCountsEndo+1)
                  LmeansEndo <- rowMeans( LCountsEndo )
                  LVarsEndo <- rowVars( LCountsEndo )
                  
                  xg = LmeansEndo
                  
                  Var_techEndo_logfit_loess =  predict(fit, LmeansEndo)
                  
                  minVar_Endo = min(LVarsEndo[LmeansEndo>2.5])
                  
                  if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)){
                    idx = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)
                    Var_techEndo_logfit_loess[idx] = 0.6*minVar_Endo       
                  }      
                  
                  is_het = (Var_techEndo_logfit_loess < LVarsEndo) &  LmeansEndo>0.3
                  res.list[["residual"]] <- LVarsEndo - Var_techEndo_logfit_loess
                  res.list[["x"]] <- LmeansEndo
                  res.list[["y"]] <- LVarsEndo
                  
                  if(plot==TRUE){

                    par(cex.lab=1.5,cex.main=1.5,mar=c(5,5,4,2))
                    plot( LmeansEndo, LVarsEndo, log="y", col=1+is_het,pch=20,cex=0.3,
                         xlab='means',ylab='Vars',main=expression("Vars ~ Means (log transformed)"))
                    xg <- seq( 0, 5.5, length.out=100 )
                    Var_techEndo_logfit_loess =  predict(fit, xg)
                    if(any(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)){
                      idx_1 = which(xg>2.5 & Var_techEndo_logfit_loess<0.6*minVar_Endo)[1]
                      idx_end = length(Var_techEndo_logfit_loess)
                      Var_techEndo_logfit_loess[idx_1:idx_end] = 0.6*minVar_Endo       
                    }
                    lines( xg, Var_techEndo_logfit_loess,lwd=2,col='green' )
                    legend('bottomright',c('Endo. genes','Var. genes',"Fit"),
                           pch=c(20,20,NA),lty = c(NA,NA,1),col=c('black','red', 'green'),cex=1.2)
                  }
                }
                is_het[is.na(is_het)] = FALSE
                res.list[["is_het"]] <- is_het
                ##res.list[["geneID"]] <- rownames(nCountsEndo)
                ###is_het
                res.list
            })


#  “memo sort” method by B. Arman Aksoy (https://gist.github.com/armish/564a65ab874a770e2c26)
#' 
#' @param M (matrix) input data
memoSort <- function(M) {
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
    scoreCol <- function(x) {
        score <- 0;
        for(i in 1:length(x)) { 
            if(x[i]) {
                score <- score + 2^(length(x)-i);
            }
        }
        return(score);
    }
    scores <- apply(M[geneOrder, ], 2, scoreCol);
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    return(M[geneOrder, sampleOrder]);
}

auto.colSet <- function(n=2,name="Set1"){
    suppressPackageStartupMessages(require("RColorBrewer"))
    if(n<=9){
         ret <- brewer.pal(max(n,3),name)[seq_len(n)]
    }else{
        ret <- colorRampPalette(brewer.pal(12,"Paired"))(n)
    }
    return(ret)
}

run.GSVA <- function(exp.dat,gset,gsva.method="gsva",bool.rnaseq=TRUE,ncores=4,min.sz=1,ssgsea.norm=T)
{
    suppressPackageStartupMessages(require("GSVA"))
    #suppressPackageStartupMessages(require("doParallel"))
    #registerDoParallel(cores = ncores)
    gsva.out <- gsva(exp.dat,gset,method=gsva.method,rnaseq=bool.rnaseq,mx.diff=T,parallel.sz=ncores,min.sz=min.sz,ssgsea.norm=ssgsea.norm)
    return(gsva.out)
}

calInfiltrationScore <- function(x,use.scale=F){
    suppressPackageStartupMessages(require("circlize"))
    suppressPackageStartupMessages(require("RColorBrewer"))
    x.scale <- t(apply(x,1,scale))
    colnames(x.scale) <- colnames(x)
    TSet <- c("T_cells","CD8_T_cells","T_helper_cells","Tcm","Tem","Th1_cells","Th2_cells","Th17_cells","TReg")
    OSet <- c(c("Macrophages","DC","B_cells","Eosinophils","Mast_cells","Neutrophils","NK_cells"),TSet)
    inf.score <- data.frame(TIS=apply(x.scale[TSet,],2,mean),
                            OIS=apply(x.scale[OSet,],2,mean),
                            CD8_Treg_Ratio=if(use.scale) apply(x.scale[c("CD8_T_cells","TReg"),],2,function(x){x[1]-x[2]}) else 
                                apply(x[c("CD8_T_cells","TReg"),],2,function(x){x[1]-x[2]}),
                            Th17_Th2_Ratio=if(use.scale) apply(x.scale[c("Th17_cells","Th2_cells"),],2,function(x){x[1]-x[2]}) else 
                                apply(x[c("Th17_cells","Th2_cells"),],2,function(x){x[1]-x[2]})
                            )
    inf.score$TIS.raw <- inf.score$TIS
    inf.score$OIS.raw <- inf.score$OIS
    q.tis <- quantile(inf.score$TIS)
    q.ois <- quantile(inf.score$OIS)
    print(q.tis)
    print(q.ois)
    #print(head(inf.score$TIS))
    #print(head(inf.score$OIS))
    for(i in 1:4){ inf.score$TIS[ q.tis[i] <= inf.score$TIS.raw & inf.score$TIS.raw <= q.tis[i+1] ] <- i }
    for(i in 1:4){ inf.score$OIS[ q.ois[i] <= inf.score$OIS.raw & inf.score$OIS.raw <= q.ois[i+1] ] <- i }
    #print(head(inf.score$TIS))
    #print(head(inf.score$OIS))
    #inf.score$TIS[inf.score$TIS > 3] <- 3
    #inf.score$TIS[inf.score$TIS < -3] <- -3
    #inf.score$OIS[inf.score$OIS > 3] <- 3
    #inf.score$OIS[inf.score$OIS < -3] <- -3
    i.col <- colorRamp2(1:4,colorRampPalette((brewer.pal(11,"BrBG")[c(7,8,9,10)]))(4), space="LAB")
    return(list(score=inf.score,
                col=list(TIS=i.col, OIS=i.col),
                x.scaled=x.scale
               ))
}

