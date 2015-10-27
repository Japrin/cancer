#!/usr/local/bin/Rscript

# excute in a dir up-level "reportData"
inFileDataProduction<-"./reportData/data.production.txt";
inPicDepth<-"./reportData/depth.png";
inFileSomaticSNVSummary<-"./reportData/somatic.snv.table";
inFileSomaticSNV<-"./reportData/all.varScan.somatic.snv.genome.final.txt";
inFileSomaticIndelSummary<-"./reportData/somatic.indel.table";
inFileSomaticIndel<-"./reportData/all.GATK.somatic.indel.genome.final.txt";
inFileSomaticCNVSummary<-"./reportData/somatic.cnv.table";
inFileSomaticCNV<-"./reportData/somatic.cnv.data";
inFileSomaticSVSummary<-"./reportData/somatic.sv.table";
inFileSomaticSV<-"./reportData/Zigongneimo.breakdancer.ann.txt";
inPicCircos<-"./reportData/circos.png";
inPicPipeline<-"./reportData/pipeline.png";
svPicDir<-"./reportData/svPic";
cnvPicDir<-"./reportData/cnvPic";
snvPicDir<-"./reportData/snvPic";
mutationSignatureDir<-"./reportData/mutationSignature";

#dir.create( "reports", showWarnings=FALSE );

require( Nozzle.R1 )

### --------------------------------- begin data -------------------------

refVarScan			<-newCitation(authors="Koboldt, D. C. etc",title="VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing",publication="Genome Research",issue="3",number="22",pages="568-76",year="2012");
refGATK				<-newCitation(authors="DePristo, M. A. etc",title="A framework for variation discovery and genotyping using next-generation DNA sequencing data",publication="Nature Genetics",issue="5",number="43",pages="491-8",year="2011");
refRSWSeq			<-newCitation(authors="Kim, T. M. etc",title="rSW-seq: algorithm for detection of copy number alterations in deep sequencing data",publication="BMC Bioinformatics",issue="",number="11",pages="432",year="2010");
refCNVSeq			<-newCitation(authors="Xie, C. etc",title="CNV-seq, a new method to detect copy number variation using high-throughput sequencing",publication="BMC Bioinformatics",issue="",number="10",pages="80",year="2009");
refBreakdancer		<-newCitation(authors="Chen, K. etc",title="BreakDancer: an algorithm for high-resolution mapping of genomic structural variation",publication="Nature Methods",issue="9",number="6",pages="677-81",year="2009");
refCircos			<-newCitation(title="Circos", url="http://circos.ca/");
### --------------------------------- end data -------------------------

# Phase 1: create report elements
r           <- newCustomReport( "癌症基因组测序分析报告" );
#r			<-setLogo(r,"./reportData/logo.png",LOGO.TOP.CENTER);

sBg         <- newSection( "背景" );
#sData       <- newSection( "实验流程" );
sData       <- newSection( "数据描述" );
sAnalyze    <- newSection( "数据分析流程" );
#sResult     <- newSection( "体细胞突变检测结果" );
sResult     <- newSection( "体细胞突变检测结果" );
sReference	<- newSection( "参考文献" );

### background ###
p           <- newParagraph( "样本编号",asStrong("NHZY0001") );
sBg         <-addTo(sBg,p);

### workflow ###
p			<-newFigure(inPicPipeline,"Fastaq格式的read数据用BWA[1]比对到参考基因组（hg19）上得到初步的BAM[2]文件；BAM文件再用Picard[3]、GATK[4]进行去重复（duplicate removal）、局部重比对（local realignment）、碱基质量值重校正（base quality recalibration）等处理，得到最终的BAM文件。\nNormal和Tumor的最终的BAM文件都用GATK的UnifiedGenotyper模块进行SNP/INDEL的检测，用CNVnator[5]进行CNV的检测，用BreakDancer[6]进行SV的检测。同时比较Normal和Tumor，可以检测体细胞突变。从VarScan[7]检测的结果过率掉群体中的多态性位点（来自dbSNP[8]）得到最后的somatic SNV；用GATK的SomaticIndelDetector模块检测somatic indel，也过滤掉群体中的多态性位点；Tumor中的拷贝数改变用rSW-Seq[9]检测；另外用breakdancer[6]检测其它体细胞结构突变。\n检测出的突变用ANNOVAR[10]进行注释；检测出的突变与肿瘤体细胞突变数据库COSMIC[11]比较；突变的基因用KEGG，GO，OMIM等数据库注释，同时与已知癌基因集比较。",fileHighRes=inPicPipeline);
sAnalyze	<-addTo(sAnalyze,p)

### data ###

table		<-read.table(inFileDataProduction,header=T,sep="\t");
p			<-newTable(table,"数据量统计");
sData		<-addTo(sData,p);
p			<-newFigure(inPicDepth,"测序深度. 深度分布（左）和累积深度分布（右）");
sData		<-addTo(sData,p);

### result ###
## somatic snv
ss    		<-newSection( "Somatic SNV" );
ss			<-addTo(ss,newParagraph("利用VarScan",asReference(refVarScan),"检测体细胞点突变（somatic SNV）。检测结果汇总如下:"));
table		<-read.table(inFileSomaticSNVSummary,header=T,sep="\t",row.names=1);
p1			<-newParagraph(paste("Exonic",asStrong(asEmph(table["Exonic",])),sep=", "))
p2			<-newParagraph(paste("splicing",asStrong(asEmph(table["splicing",])),sep=", "))
p3			<-newParagraph(paste("intronic",table["intronic",],sep=", "))
p4			<-newParagraph(paste("UTR5",table["UTR5",],sep=", "))
p5			<-newParagraph(paste("UTR3",table["UTR3",],sep=", "))
p6			<-newParagraph(paste("upstream",table["upstream",],sep=", "))
p7			<-newParagraph(paste("downstream",table["downstream",],sep=", "))
p8			<-newParagraph(paste("ncRNA",table["ncRNA",],sep=", "))
p9			<-newParagraph(paste("intergenic",table["intergenic",],sep=", "))
pp1			<-newParagraph(paste("nonsynonymous",asStrong(asEmph(table["nonsynonymousSNV",])),sep=", "))
pp2			<-newParagraph(paste("stopgain",asStrong(asEmph(table["stopgainSNV",])),sep=", "))
pp3			<-newParagraph(paste("synonymous",asStrong(asEmph(table["synonymousSNV",])),sep=", "))
pp4			<-newParagraph(paste("unknown",((table["unknown",])),sep=", "))
pp			<-newList(pp1,pp2,pp3,pp4)
p			<-newList(p1,pp,p2,p3,p4,p5,p6,p7,p8,p9)
ss			<-addTo(ss,p);


resultMutationSpectrum<-addTo(newResult("mutation spectrum",isSignificant=FALSE),
							addTo(newSection("mutation spectrum"),
								newParagraph("mutation spectrum"),
								newFigure(paste(mutationSignatureDir,"/","Spectrum.png",sep=""),"",fileHighRes=paste(mutationSignatureDir,"/","Spectrum.png",sep=""))
								)
							);
resultStrandbiasMutation<-addTo(newResult("strand bias mutation",isSignificant=FALSE),
							addTo(newSection("strandbias mutation"),
								newParagraph("strandbias mutation"),
								newFigure(paste(mutationSignatureDir,"/","Strandbias_spectrum.png",sep=""),"",fileHighRes=paste(mutationSignatureDir,"/","Strandbias_spectrum.png",sep=""))
								)
							);
resultMutationContext<-list();		
context<-c("TG","TC","TA","CA","CT","CG")
for(i in seq_along(context))  
{		
	resultMutationContext[[i]]<-addTo(newResult(context[i],isSignificant=FALSE),
							addTo(newSection("mutation context"),
								newParagraph("mutation context"),
								newFigure(paste(mutationSignatureDir,"/","mutation_context_",context[i],".png",sep=""),"",fileHighRes=paste(mutationSignatureDir,"/","mutation_context_",context[i],".png",sep=""))
								)
							);									
}
print(names(resultMutationContext))
ss			<-addTo(ss,newParagraph("Mutation signatures such as:",asSummary(resultMutationSpectrum),",",asSummary(resultStrandbiasMutation),", mutation context of ",
							asSummary(resultMutationContext[[1]]),asSummary(resultMutationContext[[2]]),asSummary(resultMutationContext[[3]]),asSummary(resultMutationContext[[4]]),asSummary(resultMutationContext[[5]]),asSummary(resultMutationContext[[6]])," etc."));

ss			<-addTo(ss,newParagraph("检测详细结果如下:"));
table		<-read.table(inFileSomaticSNV,header=T,sep="\t");
f			<-table[,"Func"]=="exonic" | table[,"Func"]=="splicing"
table		<-table[f,]
fSig		<-(grepl("nonsynonymous",table[,"ExonicFunc"]) & (table[,"AVSIFT"]<0.05 | table[,"LJB_PolyPhen2_Pred"]=="D")) | table[,"Func"]=="splicing"
fSig		<-fSig & !is.na(fSig)
table		<-table[fSig,]
p			<-newTable(table[1:min(dim(table)[1],20),],"somatic snv列表(显示了部分外显子及splicing的)",file=inFileSomaticSNV);
for(i in 1:min(dim(table)[1],20))
{	
	picFile<-dir(snvPicDir, pattern = paste(table[i,10],".*.png",sep=""), full.names=TRUE, ignore.case = TRUE)[1];
	result <- addTo( newResult( "", isSignificant=TRUE), 
					addTo( newSection("Gene Structure & MultiSequenceAlignment"), newParagraph( "Gene(protein) structure of  ", asStrong(asEmph(table[i,10]))), newFigure(picFile,"",fileHighRes=picFile) ));
	p <- addTo( p, result, row=i, column=10 );		
}
ss			<-addTo(ss,p);

ss			<-addTo(ss,newParagraph(asStrong("注:")));
pp1			<-newParagraph("CHROM=染色体");
pp2			<-newParagraph("POS=坐标");
pp3			<-newParagraph("REF=参考等位基因");
pp4			<-newParagraph("ALT=突变等位基因");
pp5			<-newParagraph("Func=exonic，intronic等");
pp6			<-newParagraph("ExonicFunc=nonsynonymous，synonymous等");
pp7			<-newParagraph("Gene=基因");
pp8			<-newParagraph("AAChange=氨基酸改变");
pp9			<-newParagraph("cosmic61=cosmic中的记录");
pp10		<-newParagraph("Conserved=是否保守区域");
pp11		<-newParagraph("SegDup=是否大片段复制区域");
ss			<-addTo(ss,newList(pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,pp10,pp11))

#AVSIFT						SIFT预测的pvaule;小于0.05为预测有害	
#LJB_SIFT					LJB中的SIFT（处理过的SIFT，与上行相反，越接近1，突变有害可能性越大
#LJB_SIFT_Pred				LJB中的SIFT预测的判断。D为有害（Damage），T为无害（Torrent）
#LJB_PolyPhen2				LJB中的PolyPhen2预测的分值
#LJB_PolyPhen2_Pred			LJB中的PolyPhen2预测的判断。D为有害（Damage），P为可能有害（Probably），B为无害（Benign）
#LJB_PhyloP					LJB中的PhyloP预测的分值
#LJB_PhyloP_Pred				LJB中的PhyloP预测的判断
#LJB_MutationTaster			LJB中的MutationTaster预测的分值
#LJB_MutationTaster_Pred		LJB中的MutationTaster预测的判断
#LJB_LRT						LJB中的LRT预测的分值
#LJB_LRT_Pred				LJB中的LRT预测的判断
#LJB_GERP++					LJB中的GERP++预测的分值
#Normal.GT					normal中的基因型
#Normal.DP					normal中的深度
#Normal.FREQ					normal中的突变等位基因的频率
#Tumor.GT					tumor中的基因型
#Tumor.DP					tumor中的深度
#Tumor.FREQ					tumor中的突变等位基因的频率
sResult		<-addTo(sResult,ss); 

## somatic indel
ss    		<-newSection( "Somatic INDEL" );
ss			<-addTo(ss,newParagraph("利用GATK",asReference(refGATK),"检测长度小于50 bp的体细胞插入缺失突变(somatic InDel)。检测结果汇总如下："));
table		<-read.table(inFileSomaticIndelSummary,header=T,sep="\t",row.names=1);
p1			<-newParagraph(paste("Exonic",asStrong(asEmph(table["Exonic",])),sep=", "))
p2			<-newParagraph(paste("splicing",asStrong(asEmph(table["splicing",])),sep=", "))
p3			<-newParagraph(paste("intronic",table["intronic",],sep=", "))
p4			<-newParagraph(paste("UTR5",table["UTR5",],sep=", "))
p5			<-newParagraph(paste("UTR3",table["UTR3",],sep=", "))
p6			<-newParagraph(paste("upstream",table["upstream",],sep=", "))
p7			<-newParagraph(paste("downstream",table["downstream",],sep=", "))
p8			<-newParagraph(paste("ncRNA",table["ncRNA",],sep=", "))
p9			<-newParagraph(paste("intergenic",table["intergenic",],sep=", "))
pp1			<-newParagraph(paste("Nonframeshift deletion",asStrong(asEmph(table["Nonframeshift deletion",])),sep=", "))
pp2			<-newParagraph(paste("Nonframeshift insertion",asStrong(asEmph(table["Nonframeshift insertion",])),sep=", "))
pp3			<-newParagraph(paste("frameshift deletion",asStrong(asEmph(table["frameshift deletion",])),sep=", "))
pp4			<-newParagraph(paste("frameshift insertion",asStrong(asEmph(table["frameshift insertion",])),sep=", "))
pp5			<-newParagraph(paste("unknown",((table["unknown",])),sep=", "))
pp			<-newList(pp1,pp2,pp3,pp4,pp5)
p			<-newList(p1,pp,p2,p3,p4,p5,p6,p7,p8,p9)
ss			<-addTo(ss,p);
ss			<-addTo(ss,newParagraph("检测详细结果如下:"));
table		<-read.table(inFileSomaticIndel,header=T,sep="\t");
f			<-table[,"Func"]=="exonic" | table[,"Func"]=="splicing"
table		<-table[f,]
p			<-newTable(table[c(1:min(dim(table)[1],20)),],"somatic indel列表(显示了外显子及splicing的)",file=inFileSomaticIndel);
for(i in 1:min(dim(table)[1],20))
{	
	result <- addTo( newResult( "", isSignificant=TRUE), 
					addTo( newSection("Gene Structure & MultiSequenceAlignment"), newParagraph( "Gene(protein) structure of  ", asStrong(asEmph(table[i,10]))) ) );
	p <- addTo( p, result, row=i, column=10 );		
}
ss			<-addTo(ss,p);
sResult		<-addTo(sResult,ss);
### somatic cnv
ss    		<-newSection( "Somatic CNV" );
ss			<-addTo(ss,newParagraph("利用CNVseq",asReference(refCNVSeq),"检测拷贝数改变（somatic CNV）。检测结果汇总如下："));
table		<-read.table(inFileSomaticCNVSummary,header=T,sep="\t",row.names=1);
p1			<-newParagraph(paste("total",table["total",],sep=", "))
pp1			<-newParagraph(paste("gain",table["gain",],sep=", "))
pp2			<-newParagraph(paste("loss",table["loss",],sep=", "))
ss			<-addTo(ss,newList(p1,newList(pp1,pp2)))
p1			<-newParagraph(paste("total affected genes",table["cna_gene",],sep=", "))
pp1			<-newParagraph(paste("\"gain\" gene",table["gain_gene",],sep=", "))
pp2			<-newParagraph(paste("\"loss\" gene",table["loss_gene",],sep=", "))
ss			<-addTo(ss,newList(p1,newList(pp1,pp2)))
ss			<-addTo(ss,newParagraph("检测详细结果如下:"));
table		<-read.table(inFileSomaticCNV,header=T,sep="\t");
p			<-newTable(table[c(1:10),],"somatic CNV列表(显示了前十个)",file=inFileSomaticCNV);
for(i in 1:min(dim(table)[1],10))
{
	chr<-table[i,1]
	id<-table[i,7]
	picFile<-paste(cnvPicDir,"/",chr,"/",chr,"_",id,".png",sep="")
	result <- addTo( newResult( "", isSignificant=FALSE), 
					addTo( newSection("Gene Structure & CNV"), newParagraph( "Detail of CNV  ", asStrong(asEmph(table[i,7]))), newFigure(picFile,"\"refGene\"显示的是基因的结构，\"RegionDetail\"中的淡绿色块指示cnv的区域，细节显示中了每个窗口的log2 ratio。",fileHighRes=picFile) ) );
	p <- addTo( p, result, row=i, column=7 );
}
ss			<-addTo(ss,p);
sResult		<-addTo(sResult,ss); 

### somatic sv
ss    		<-newSection( "Somatic SV" );
table		<-read.table(inFileSomaticSVSummary,header=T,sep="\t");
ss			<-addTo(ss,newParagraph("利用breakdancer",asReference(refBreakdancer),"检测大的体细胞结构突变（somatic SV）。全基因组范围内共检测出",asStrong(asEmph(sum(table[,-1]))),"个somatic SV"));
p			<-newTable(table,"somatic SV汇总");
ss			<-addTo(ss,p);
table		<-read.table(inFileSomaticSV,header=T,sep="\t",stringsAsFactors=FALSE);
f			<-!grepl("intergenic",table[,"region"])
table		<-table[f,]
p			<-newTable(table[c(1:10),],"somatic SV列表(显示了前十个)",file=inFileSomaticSV);
for(i in 1:min(dim(table)[1],10))
{	
	
	m<-regexec("SVID ([0-9]+)", table[i,8])
	sv_id<-(regmatches(table[i,8], m)[[1]][2])
	pattern<-paste("[.]",sv_id,"[.]png",sep="")
	picFile<-dir(svPicDir, pattern = pattern, full.names=TRUE, ignore.case = TRUE)[1]
	result <- addTo( newResult( "", isSignificant=FALSE), 
					addTo( newSection("Gene Structure & SV"), newParagraph( "Details: sv affect SV ", asStrong(asEmph(table[i,4]))), newFigure(picFile,"\"refGene\"显示的是基因的结构，\"RegionDetail\"中的淡绿色块指示的refGene轨道中将显示细节的区域，细节显示中，read用淡蓝色的块表示。",fileHighRes=picFile) ) );
	p <- addTo( p, result, row=i, column=8 );		
}
ss			<-addTo(ss,p);
sResult		<-addTo(sResult,ss);

### circos
ss    		<-newSection( "circos" );
ss			<-addTo(ss,newParagraph("所有的体细胞突变用Circos",asReference(refCircos),"图展示"));
ss			<-addTo(ss,newFigure(inPicCircos,"用Circos展示所有的体细胞突变. 从外到内的轨道依次为染色体，基因密度，拷贝数变化（红色为loss，蓝色为gain），插入，倒位（橙色), 缺失（蓝色），10M窗口的somatic snv数，10M窗口的somatic indel数，染色体内重排（绿色），染色体间重排（紫色）"));
sResult		<-addTo(sResult,ss);
# reference

sReference	<-addTo(sReference, refVarScan,refGATK,refCNVSeq,refBreakdancer,refCircos);


# Phase 2: assemble report structure bottom-up
r <- addTo( r, sBg, sAnalyze, sData, sResult, sReference);

# Phase 3: render report to file
writeReport( r, filename="my_report" ); # w/o extension
