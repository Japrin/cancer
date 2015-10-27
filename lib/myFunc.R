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

qcAndVizMat<-function(vstMat,designM,outDir,colSet=NULL,intgroup=c("sampleType"),ntop=500,extra="",sfilter=NULL, gfilter=NULL,complexHeatmap.use=NULL,...)
{
    if(!is.null(sfilter))
    {
	    vstMat <- vstMat[,sfilter, drop = FALSE]
	    designM <- designM[sfilter,, drop = FALSE]
    }
    if(!is.null(gfilter))
    {
	    vstMat <- vstMat[rownames(vstMat) %in% gfilter[!is.na(gfilter)],, drop = FALSE]
    }
	vstMat <- vstMat[, rownames(designM), drop = FALSE]

    if(is.null(colSet))
    {
        nLevel <- length(levels(designM$sampleType))
        if(nLevel<=9)
        {
            suppressPackageStartupMessages(require("RColorBrewer"))
            colSet <- brewer.pal(nLevel,"Set1")
            names(colSet) <- levels(designM$sampleType)
            patientcolors <- colSet[as.numeric(designM$sampleType)]
        }else
        {
            colSet <- rainbow(nLevel)
            names(colSet) <- levels(designM$sampleType)
            patientcolors <- colSet[as.numeric(designM$sampleType)]
        }
    }else
    {
      colSet <- colSet[names(colSet) %in% levels(designM$sampleType)]
      patientcolors <- colSet[as.character(designM$sampleType)]
    }
    ## PCA
    #require("DESeq2")
    suppressPackageStartupMessages(require("ggplot2"))
    suppressPackageStartupMessages(require("gplots"))
    suppressPackageStartupMessages(require("factoextra"))
    rowVar <- apply(vstMat,1,var)
    #rowVar <- apply(vstMat,1,function(x){ var(x)/(mean(x)^2) } )
    if(is.null(ntop))
    {
        ntop=length(rowVar)
    }
    if(ntop<3)
    {
        loginfo(paste0("Too few samples: n=",ntop))
        return
    }
    cat(sprintf("%s\t---- PCA plot\n", Sys.time()))
    select <- order(rowVar, decreasing = TRUE)[seq_len(min(ntop, length(rowVar)))]
    dat.plot <<- t(vstMat[select, ])
    cnames <- entrezToXXX(colnames(dat.plot))
    cnames.na <- which(is.na(cnames))
    cnames[cnames.na] <- colnames(dat.plot)[cnames.na]
    colnames(dat.plot) <- cnames

    pca <<- prcomp(dat.plot,scale=T)
    percentVar <<- pca$sdev^2/sum(pca$sdev^2)

    pca.eig <- (pca$sdev)^2
    variance.in.percent <- pca.eig*100/sum(pca.eig)
    cumvar.in.percent <- cumsum(variance.in.percent)
    pca.eig.df <- data.frame(eig = pca.eig, variance = variance.in.percent, cumvariance = cumvar.in.percent)
    #eig.val <- get_eigenvalue(pca)
    pca.var <- get_pca_var(pca)
    pca.ind <- get_pca_ind(pca)

    d <- data.frame(geneSymbol=rownames(pca.var$contrib))
    d <- cbind(d,pca.var$contrib)
	write.table(d, paste(outDir,"/PCA",extra,".contribution.txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = F)

    if(!all(intgroup %in% names(designM))) {
        stop("the argument 'intgroup' should specify columns of designM")
    }
    intgroup.df <- as.data.frame(designM[, intgroup, drop = FALSE])
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, intgroup.df, names = colnames(vstMat))
	write.table(d, paste(outDir,"/PCA",extra,".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = F)
    
    print(head(percentVar)) 
	
    pdf(file=paste(outDir,"/PCA",extra,".pdf",sep=""),width=10,height=8)
    par(mar=c(5,5,4,8),cex.lab=1.5)
    ### scree plot
    nDim.plot <- min(30,nrow(pca.eig.df))
    xat <- barplot(pca.eig.df[1:nDim.plot, "variance"], names.arg=1:nDim.plot, main = "Variances", xlab = "Principal Components", ylab = "Percentage of variances", col ="steelblue")
    lines(x = xat, pca.eig.df[1:nDim.plot, "variance"], type="b", pch=19, col = "red")
    abline(h = 1*100/sum(pca.eig), lty = 2, col = "red", lwd=2)

    ### loading plot
    # Plot the correlation circle
    #a <- seq(0, 2*pi, length = 100)
    #plot( cos(a), sin(a), type = 'l', col="gray", xlab = "PC1",  ylab = "PC2")
    #abline(h = 0, v = 0, lty = 2)
    #arrows(0, 0, pca.var$coord[, 1], pca.var$coord[, 2], length = 0.1, angle = 15, code = 2)
    #text(pca.var$coord[, 1], pca.var$coord[, 2], labels=rownames(pca.var$coord), cex = 1, adj=1)
    ### score plot
    plot(x=NULL,y=NULL, xlim = range(d$PC1), ylim = range(d$PC2), type = "n", main=sub("^\\.","",extra,perl=T),cex.main=1.5,
         xlab=paste0("PC1: ", round(percentVar[1] * 100,2), "% variance"), 
         ylab=paste0("PC2: ", round(percentVar[2] * 100,2), "% variance"))
    invisible(lapply(levels(d$group),function(x){ points(subset(d,group==x,select="PC1")[,1],subset(d,group==x,select="PC2")[,1],col=colSet[x],pch=16,cex=1.2)  }))
    legend("right",legend=names(colSet),fill = NULL,inset = -0.19,xpd = NA,cex=1.5,pch=16,border =NA,col = colSet)

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
    print(fviz_pca_ind(pca, label = "none", col.ind = patientcolors ,habillage = group, addEllipses = T) + theme_minimal())

    dev.off()


    # NMF
    dat.plot <- t(dat.plot)
    doNMF <- function(dat.plot)
    {
        suppressPackageStartupMessages(require("NMF"))
        nmf.res <<- nmf(dat.plot[,], 2:7, nrun = 100, .opt = "vp8", seed = 123456)
        pdf(file=paste(outDir,"/NMF",extra,".pdf",sep=""),width=12,height=10)
        print(plot(nmf.res))
        ###nmf.options(grid.patch=TRUE)
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         consensusmap(x[[i]],
                                      annCol = designM[,"sampleType",drop=F],
                                      annColors = list(sampleType=colSet),
                                      fontsize=20) 
                    }, nmf.res$fit))
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         basismap(x[[i]],main = paste0("Metagenes (Rank=",i+1,")"),fontsize=20) 
                    }, nmf.res$fit))
        invisible(sapply(seq(length(nmf.res$fit)), function(i,x){ 
                         coefmap(x[[i]],main = paste0("Metagene contributions in each sample (Rank=",i+1,")"),
                                 annCol = designM[,"sampleType",drop=F],
                                 annColors = list(sampleType=colSet),
                                 fontsize=20) 
                    }, nmf.res$fit))
        dev.off()
        save(nmf.res,file = paste(outDir,"/NMF",extra,".RData",sep=""))
    }
    tryCatch(doNMF(dat.plot), error = function(e) e, finally = loginfo(sprintf("NMF run finished")) )
    ## t-SNE
    library(Rtsne)
    dat.plot <- t(vstMat[select, ])
    doit <- function(par.perplexity)
    {
        Rtsne.res <- Rtsne(dat.plot[,],perplexity=par.perplexity)
        pdf(file=sprintf("%s/tsne%s.perplexity%d.pdf",outDir,extra,par.perplexity),width=10,height=8)
        par(mar=c(5,5,4,8),cex.lab=1.5,cex.main=1.5)
        plot(Rtsne.res$Y, t='n', main="BarnesHutSNE",xlab="Dim1",ylab="Dim2")
        #text(Rtsne.res$Y, labels=rownames(dat.plot)) 
        points(Rtsne.res$Y,col=colSet[as.character(designM$sampleType)],pch=16)
        legend("right",legend=names(colSet),fill = NULL,inset = -0.19,xpd = NA,cex=1.5,pch=16,border =NA,col = colSet)
        dev.off()
    }
    sapply(seq(5,50,5), function(i) { tryCatch(doit(i), error = function(e) e, finally = loginfo(sprintf("Rtsne run finished with perplexity %d \n",i)) ) } )


    #print(patientcolors)
    designM$color <- patientcolors
    legendDat <- unique(designM[,c("sampleType","color")])
    
    #for(n in c(50,100,150,200,250,300,350,400,450,500,1000,2000,3000,4000,5000))
    for(n in c(50,100,150,200,250,300,350,400,450,500,1000,2000))
    {
        ###hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
        select <- order(rowVar,decreasing = TRUE)[seq_len(min(n, length(rowVar)))]
        dat.plot <-as.matrix(vstMat[select,])
        
        dat.plot.export <- data.frame(geneID=rownames(dat.plot),geneSymbol=entrezToXXX(rownames(dat.plot)))
        dat.plot.export <- cbind(dat.plot.export,dat.plot)
	    write.table(dat.plot.export,paste(outDir,"/QC.Cluster.Var.n",n,extra,".txt",sep=""),sep="\t",col.names = T,row.names = F,quote = F)
        
        nn <- min(n, length(rowVar))
        m <- ncol(dat.plot)
        rownames(dat.plot) <- entrezToXXX(rownames(dat.plot))
        if(is.null(complexHeatmap.use))
        {
            pdf(paste(outDir,"/QC.Cluster.Var.n",n,extra,".pdf",sep=""),width=12,height=12)
            ###heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20), main="Most variable genes",cexRow=min(1.8,55/nn),cexCol=min(1.8,55/m),distfun=function(x){ as.dist(1-cor(t(x),method = "spearman")) } )
            heatmap.2(dat.plot,col=bluered(100), ColSideColors=patientcolors, Rowv = T, Colv = T, scale="row", density.info="none", dendrogram="both", keysize=1.2, trace="none", margin=c(15, 20), main="Most variable genes",cexRow=min(1.8,55/nn),cexCol=min(1.8,55/m))
            legend("topright",legend=legendDat$sampleType,fill=legendDat$color,border=legendDat$color,cex=1.5)
            dev.off()
        }else
        {
            require("ComplexHeatmap")
            require("circlize")
            pdf(paste(outDir,"/QC.Cluster.Var.n",n,extra,".pdf",sep=""),width=15,height=12)

            annDF <- data.frame(sampleType=designM$sampleType)
            annColList <- list(sampleType=colSet)
            ha.col <- HeatmapAnnotation(df = annDF, col = annColList, show_legend = TRUE )
            top_annotation_height <- unit(1.5, "cm")
            ht <- Heatmap(dat.plot, name="NormalizedRC",
                      col = colorRamp2(c(0, 0.609, 1, 10), c("darkblue", "darkblue", "yellow", "red")),
                      column_hclust_height = unit(6, "cm"), row_hclust_width = unit(6, "cm"),
                      column_names_gp = gpar(fontsize = 12*55/m),row_names_gp = gpar(fontsize = 12*55/max(n,25)),
                      show_heatmap_legend = T, row_names_max_width = unit(10,"cm"),
                      top_annotation_height = top_annotation_height,
                      clustering_distance_rows = "spearman", clustering_distance_columns = "spearman",  ### euclidean
                      top_annotation = ha.col,...)
            draw(ht, legend_grid_width = unit(0.8, "cm"), legend_grid_height = unit(0.8, "cm"),
                legend_title_gp = gpar(fontsize = 14, fontface = "bold"),
                legend_label_gp = gpar(fontsize = 12))
            dev.off()
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
    cat(sprintf("%s\t---- qcAndViz finished\n", Sys.time()))

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
    in.table <- read.table(in.file,header = T,check.names = F,stringsAsFactors = F,sep="\t",comment.char="!")
    ret <- in.table$color
    names(ret) <- in.table$sampleType
    return(ret)
}


