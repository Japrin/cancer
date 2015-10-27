#!/PROJ/GR/share/Software/R/bin/Rscript

args<-commandArgs(T)
if(length(args)<2)
{
	cat("usage: plot.breakdancer.R <infile> <prefix>\n")
	q()
}
library("grid")

infile<-args[1]
prefix<-args[2]
a<-read.table(infile,sep="\t",header=F)
colnames(a)<-c("rid","chr1","pos1","chr2","pos2","isize")
#HWI-ST1276:121:C28LNACXX:5:1306:3666:5108       1       112691564       1       112704731       13267

CHR<-a$chr1[1]
#png(paste(prefix,".png",sep=""),width=800,height=600)

X <- matrix(rexp(2000), ncol=2)

aa<-1:100
bb<-rep(50,100)
aa<-a$pos1
bb<-rep(50,length(aa))

pushViewport(plotViewport(c(5.1, 4.1, 4.1, 2.1)))
pushViewport(dataViewport(aa,bb,yscale=c(0,100)))
#grid.points(X[,1], X[,2], pch=13) # points
grid.rect() # bounding rectangle
grid.xaxis() # x-axis
grid.yaxis() # y-axis
grid.points(aa,bb,pch=16,gp=gpar(col="red"))
#grid.text(expression(italic(X[1])), y=unit(-3, "lines")) # x-axis label
#grid.text(expression(italic(X[2])), x=unit(-3, "lines"), rot=90) # y-axis label
#grid.text("Plot 1", x=0.86, y=0.9, gp=gpar(fontface="bold", cex=1.6)) # add label
#upViewport
for(i in 1:length(a$rid))
{
	x<-c(a[i,"pos1"],a[i,"pos1"],a[i,"pos2"],a[i,"pos2"])
	y<-c(25,30,30,25)
	pushViewport(dataViewport(x,y))
	grid.bezier(x,y)
	upViewport()
	#print(x)
	#vp<-viewport(x=unit((a[i,"pos1"]+a[i,"pos2"])/2,"points"),y=unit(25,"points"),width=unit(abs(a[i,"pos2"]-a[i,"pos1"]),"points"),height=unit(5,"points"))
	#print(x=unit((a[i,"pos1"]+a[i,"pos2"])/2,"points"))
	#pushViewport(vp)
	#x<-c(0,0,1,1)
	#y<-c(0,1,1,0)
	#print(x)
	#grid.bezier(x,y)
	#grid.rect()
	#upViewport()

}
upViewport()
dev.off()
