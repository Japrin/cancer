#!/usr/bin/env Rscript

#library("asbio")

#in.file <- "/WPS1/zhenglt/work/MHC/TCGA/OUT.netMHC-3.4.A02/HLA-A02.IC50"
#out.file <- "/WPS1/zhenglt/work/MHC/TCGA/OUT.netMHC-3.4.A02/HLA-A02.IC50.correlation.png"
in.file <- "/WPS1/zhenglt/work/MHC/TCGA/OUT.netMHC-3.4.A02/all.IC50"
out.prefix <- "/WPS1/zhenglt/work/MHC/TCGA/OUT.netMHC-3.4.A02/all.IC50.correlation"

data <- read.table(in.file,header = T,row.names = 1,check.names = F)
data <- data[,grep("IC50",names(data),perl = T,value = T)]
names(data) <- gsub("_IC50","",gsub("HLA-","",names(data)))

head(data)
dim(data)

source("/Share/BP/zhenglt/02.pipeline/cancer/lib/myFunc.R")

my.func<-function (x, y, col ='blue', bg = NA, pch = par("pch"), cex = 0.01, col.line = "blue", lty = par("lty"),pty=par("pty"),...) 
{
	#points(x, y, pch = ".", col = col, bg = bg, cex = cex)

	upar <- par(usr = c(0, max(log(x)+3), 0, max(log(y)+3) ))
	on.exit(par(upar))
	
	smoothScatter(log(x), log(y), add = TRUE)
	#smoothScatter(x, y, add = TRUE)
	f <- is.finite(y) & is.finite(x)
	abline(lm(y~x,data=data.frame(x=log(x[f]),y=log(y[f]))))
	#axis(3)
	#abline(a=2,b=1,col="red", bg = bg, cex = cex)
	#abline(a=-2,b=1,col="red", bg = bg, cex = cex)
	#abline(a=0,b=1,col=col, bg = bg, cex = cex)
}

plotOne <- function(dat,png.file)
{
	png(png.file,width = 800,height = 800)
	pairs2(dat,cex.labels=2,cex=2,gap=0.8,lower.panel=panel.cor,upper.panel=my.func,cex.axis=1.5,xlim=c(0,15),ylim=c(0,15))
	dev.off()
}

grep("A02:",names(data),perl = T,value = T)
plotOne(data[,grep("A02:",names(data),perl = T,value = T)],paste0(out.prefix,".A02.png"))
grep("B15:",names(data),perl = T,value = T)
plotOne(data[,grep("B15:",names(data),perl = T,value = T)],paste0(out.prefix,".B15.png"))
grep("A26:",names(data),perl = T,value = T)
plotOne(data[,grep("A26:",names(data),perl = T,value = T)],paste0(out.prefix,".A26.png"))
grep("A32:",names(data),perl = T,value = T)
plotOne(data[,grep("A32:",names(data),perl = T,value = T)],paste0(out.prefix,".A32.png"))
grep("A68:",names(data),perl = T,value = T)
plotOne(data[,grep("A68:",names(data),perl = T,value = T)],paste0(out.prefix,".A68.png"))
grep("B08:",names(data),perl = T,value = T)
plotOne(data[,grep("B08:",names(data),perl = T,value = T)],paste0(out.prefix,".B08.png"))
grep("B40:",names(data),perl = T,value = T)
plotOne(data[,grep("B40:",names(data),perl = T,value = T)],paste0(out.prefix,".B40.png"))
grep("A24:",names(data),perl = T,value = T)
plotOne(data[,grep("A24:",names(data),perl = T,value = T)],paste0(out.prefix,".A24.png"))
grep("A30:",names(data),perl = T,value = T)
plotOne(data[,grep("A30:",names(data),perl = T,value = T)],paste0(out.prefix,".A30.png"))
grep("B27:",names(data),perl = T,value = T)
plotOne(data[,grep("B27:",names(data),perl = T,value = T)],paste0(out.prefix,".B27.png"))
grep("B35:",names(data),perl = T,value = T)
plotOne(data[,grep("B35:",names(data),perl = T,value = T)],paste0(out.prefix,".B35.png"))
grep("B44:",names(data),perl = T,value = T)
plotOne(data[,grep("B44:",names(data),perl = T,value = T)],paste0(out.prefix,".B44.png"))
grep("C07:",names(data),perl = T,value = T)
plotOne(data[,grep("C07:",names(data),perl = T,value = T)],paste0(out.prefix,".C07.png"))
