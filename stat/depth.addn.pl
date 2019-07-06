#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

=pod

=head1 Usage
	perl $0 [option] <bam> <outdir>
		-q	base quality [default 0]
		-Q	mapping quality [default 0]
		-n	addn file [default "/PUBLIC/database/HEALTH/genome/human/b37_gatk/b37.84.NBlock"; "-" if don't account for N region]
		-h	help

=cut
		

my ($basethres,$mapQthres,$total_chr,$addnFile,$help);
GetOptions("q:i"=>\$basethres,"Q:i"=>\$mapQthres,"l:i"=>\$total_chr,"n:s"=>\$addnFile,"h"=>\$help);

#$total_chr ||= 2471805657;
$basethres ||= 0;
$mapQthres ||= 0;
if(!defined($addnFile) || $addnFile ne "-") { $addnFile="/WPSnew/zhenglt/02.pipeline/cancer/stat/b37.NBlock.larger20bp.bed"; }
die `pod2text $0` if(@ARGV<2 || $help);

my $bam=shift;
my $outdir=shift;

##calcualte genome length
$total_chr=`samtools view -H $bam | awk '/LN:/' | sed 's/:/\t/g' | awk '{sum+=\$NF} END{print sum}'`;
chomp $total_chr;
if($addnFile ne "-")
{
	my $_ttt=`awk '{sum+=\$4}END{print sum}' $addnFile`;
	chomp $_ttt;
	$total_chr-=$_ttt;
}
#

`mkdir -p $outdir` unless -d $outdir;

open DEPTH,"samtools depth -q $basethres -Q $mapQthres $bam | " or die;

my %depth=();
my $maxCov=0;
my $Average_sequencing_depth=0;
my $Average_sequencing_depth4=0;
my $Average_sequencing_depth10=0;
my $Average_sequencing_depth20=0;
my $Coverage=0;
my $Coverage4=0;
my $Coverage10=0;
my $Coverage20=0;

my $Coverage_bases=0;
my $Coverage_bases_4=0;
my $Coverage_bases_10=0;
my $Coverage_bases_20=0;

my $total_Coverage_bases=0;
my $total_Coverage_bases_4=0;
my $total_Coverage_bases_10=0;
my $total_Coverage_bases_20=0;

while(<DEPTH>)
{
	chomp;
	my @arr = split;
	$depth{$arr[2]}+=1;
}
close DEPTH;

my @depth=sort {$a<=>$b} keys %depth;

open HIS,">$outdir/depth_frequency.xls" or die;
open CUM,">$outdir/cumu.xls" or die;
print CUM "Depth\tPercent\n";
foreach my $depth1 (sort {$a<=>$b} keys %depth)
{
	next if($depth1==0);
	my $per=$depth{$depth1}/$total_chr;
	$total_Coverage_bases += $depth1*$depth{$depth1};
	$Coverage_bases += $depth{$depth1};

	if($depth1>=4)	
	{
		$total_Coverage_bases_4 += $depth1 * $depth{$depth1};
		$Coverage_bases_4 += $depth{$depth1};
	}
	if($depth1>=10)
	{
		$total_Coverage_bases_10 += $depth1 * $depth{$depth1};
		$Coverage_bases_10 += $depth{$depth1};
	}
	if($depth1>=20)
	{
		$total_Coverage_bases_20 += $depth1 * $depth{$depth1};
		$Coverage_bases_20 += $depth{$depth1};
	}



	$maxCov=$per if($maxCov<$per);
	my $tmp=0;
	print HIS "$depth1\t$per\t$depth{$depth1}\n";
	foreach my $depth2(@depth)
	{
		$tmp+=$depth{$depth2} if($depth2 >= $depth1); 
	}
	$tmp=$tmp/$total_chr;
	print CUM "$depth1\t$tmp\n";
}

$Average_sequencing_depth=$total_Coverage_bases/$total_chr;
$Coverage=$Coverage_bases/$total_chr;
$Average_sequencing_depth4=$total_Coverage_bases_4/$total_chr;
$Coverage4=$Coverage_bases_4/$total_chr;
$Average_sequencing_depth10=$total_Coverage_bases_10/$total_chr;
$Coverage10=$Coverage_bases_10/$total_chr;
$Average_sequencing_depth20=$total_Coverage_bases_20/$total_chr;
$Coverage20=$Coverage_bases_20/$total_chr;


print "Average_sequencing_depth:\t",sprintf("%.2f",$Average_sequencing_depth),"\n";
print "Coverage:\t",sprintf("%.2f%%",100*$Coverage),"\n";
print "Coverage_at_least_4X:\t",sprintf("%.2f%%",100*$Coverage4),"\n";
print "Coverage_at_least_10X:\t",sprintf("%.2f%%",100*$Coverage10),"\n";
print "Coverage_at_least_20X:\t",sprintf("%.2f%%",100*$Coverage20),"\n";


close HIS;
close CUM;

if(1)
{
	my $ylim = 100*$maxCov;
	my ($xbin,$ybin);
	$ylim= int($ylim) + 1;
	if($ylim <= 3)
	{
		$ybin = 0.5;
	}else{
		$ybin=1;
	}
	my $xlim=0;
	if($Average_sequencing_depth<30)
	{
		$xlim=100;
		$xbin=20;
	}elsif($Average_sequencing_depth < 50)
	{
		$xlim=160;
		$xbin=20;
	}elsif($Average_sequencing_depth  < 120)
	{
		$xlim=250;
		$xbin=50;
	}else{
		$xlim=600;
		$xbin=100;
	}
	histPlot($outdir,"$outdir/depth_frequency.txt",$ylim,$ybin,$xlim,$xbin);
	cumuPlot($outdir,"$outdir/cumu.txt",$xlim,$xbin);
}

sub cumuPlot {
	my ($outdir, $dataFile, $xlim, $xbin) = @_;
	my $figFile = "$outdir/cumuPlot.pdf";
	my $Rline=<<Rline;
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	opar <- par()
	x <- rt\$V1[1:($xlim+1)]
	y <- 100*rt\$V2[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/cumuPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,100,by=20)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
	
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("R CMD BATCH  $figFile.R");
}


sub histPlot {
	my ($outdir, $dataFile, $ylim, $ybin, $xlim, $xbin) = @_;
	my $figFile = "$outdir/histPlot.pdf";
	my $Rline=<<Rline;
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	opar <- par()
	t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
	y=c(rt\$V2[1:$xlim],t)
	y <- y*100
	x <- rt\$V1[1:($xlim+1)]
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/histPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="blue",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("R CMD BATCH  $figFile.R");
#	system("rm  $figFile.R  $figFile.Rout");
}
