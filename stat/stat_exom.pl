#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my ($bamfile,$outdir,$regionfile,$Plot,$bin,$help);

GetOptions(
		"i:s"=>\$bamfile,
		"o:s"=>\$outdir,
		"r:s"=>\$regionfile,
		"b:s"=>\$bin,
		"h"=>\$help,
		"plot"=>\$Plot,
		);

my $usage=<<USAGE;
usage:perl $0
		-i <bamfile>
		-r <region file>
		-o <outdir>
		-plot 
		-h help
USAGE

die $usage if (!$bamfile || $help || !$outdir || !$regionfile);

my $Initial_bases_on_target=0;
my $Initial_bases_near_target=0;
my $Initial_bases_on_or_near_target=0;
my $Total_effective_reads=0;
my $Total_effective_yield=0;
my $Average_read_length=0;
my $Effective_sequences_on_target=0;
my $Effective_sequences_near_target=0;
my $Effective_sequences_on_or_near_target=0;
my $Fraction_of_effective_bases_on_target=0;
my $Fraction_of_effective_bases_on_or_near_target=0;
my $Average_sequencing_depth_on_target=0;
my $Average_sequencing_depth_near_target=0;
my $Mismatch_rate_in_target_region=0;
my $Mismatch_rate_in_all_effective_sequence=0;
my $Base_covered_on_target=0;
my $Coverage_of_target_region=0;
my $Base_covered_near_target=0;
my $Coverage_of_flanking_region=0;
my $Fraction_of_target_covered_with_at_least_20x=0;
my $Fraction_of_target_covered_with_at_least_10x=0;
my $Fraction_of_target_covered_with_at_least_4x=0;
my $Fraction_of_flanking_region_covered_with_at_least_20x=0;
my $Fraction_of_flanking_region_covered_with_at_least_10x=0;
my $Fraction_of_flanking_region_covered_with_at_least_4x=0;


open BAM,"samtools view -F 0x0004 -X $bamfile | " or die $!;
open REG,"$regionfile" or die $!;
`mkdir -p $outdir` unless -d $outdir;
my %hash=();
open FREG,">$outdir/freg_tmp.txt" or die $!;
while(<REG>)
{
	chomp;
	next if(/^$/);
	my @info=split (/\t/,$_);
	$hash{$info[0]}+=$info[2]-$info[1];
	$Initial_bases_on_target+=$info[2]-$info[1];
	my $pri_beg_pos=$info[1]-200;
	my $pri_end_pos=$info[1];
	my $next_beg_pos=$info[2];
	my $next_end_pos=$info[2]+200;
	if($pri_beg_pos<0) { $pri_beg_pos=0; }
	#if ($info[0] eq 'chrM')
	#{
	#	print FREG "$info[0]\t$next_beg_pos\t$next_end_pos\n";
	#	next;
	#}
	print FREG "$info[0]\t$pri_beg_pos\t$pri_end_pos\n";
	print FREG "$info[0]\t$next_beg_pos\t$next_end_pos\n";
}
close(FREG);
close(REG);

`sortBed -i $outdir/freg_tmp.txt | mergeBed -i stdin | subtractBed -a stdin -b $regionfile   >$outdir/freg.txt`;
#`rm $outdir/freg_tmp.txt`;

$Initial_bases_near_target=`awk '{total+=\$3-\$2+1};END{print total}' $outdir/freg.txt`;
chomp($Initial_bases_near_target);

#`samtools mpileup -q 0 -Q 0 -A -B -d 1000000 --ff 0x400 -l $regionfile $bamfile | cut -f 1,2,4 > $outdir/target.depth`;
#`samtools mpileup -q 0 -Q 0 -A -B -d 1000000 --ff 0x400 -l $outdir/freg.txt $bamfile | cut -f 1,2,4 > $outdir/flanking.depth`;
`samtools depth -b $regionfile $bamfile  >$outdir/target.depth`;
`samtools depth -b $outdir/freg.txt $bamfile  >$outdir/flanking.depth`;

$Effective_sequences_on_target=`awk '{total+=\$3};END{print total}' $outdir/target.depth`;
chomp($Effective_sequences_on_target);
$Effective_sequences_near_target=`awk '{total+=\$3};END{print total}' $outdir/flanking.depth`;
chomp($Effective_sequences_near_target);
$Effective_sequences_on_or_near_target=$Effective_sequences_on_target+$Effective_sequences_near_target;

$Initial_bases_on_or_near_target=$Initial_bases_on_target+$Initial_bases_near_target;
	
my $Mismatch_base_in_target_region=0;
my $Mismatch_base_in_all_effective_sequence=0;
while(<BAM>)
{
	chomp;
	my @_F=split /\t/;
	if($_F[1]=~/d/) { next; }
	$Total_effective_reads++;
	if($_=~/XM:i:(\d+)/)
	{
		$Mismatch_base_in_all_effective_sequence+=$1;
	}
}

#`samtools mpileup -q 0 -Q 0 -A -B -d 1000000 --ff 0x400 $bamfile | cut -f 1,2,4 > $outdir/whole_genome.depth`;
`samtools depth $bamfile >$outdir/whole_genome.depth`;
$Total_effective_yield=`awk '{total+=\$3};END{print total}' $outdir/whole_genome.depth`;
chomp $Total_effective_yield;
$Average_read_length=$Total_effective_yield/$Total_effective_reads;

$Fraction_of_effective_bases_on_target=$Effective_sequences_on_target/$Total_effective_yield;
$Fraction_of_effective_bases_on_or_near_target=$Effective_sequences_on_or_near_target/$Total_effective_yield;
	
open TMP,"samtools view -F 0x0004 -X -L $outdir/freg.txt $bamfile | " or die $!;
while(<TMP>)
{
	chomp;
	my @_F=split /\t/;
	if($_F[1]=~/d/) { next; }
	if($_=~/XM:i:(\d+)/)
	{
		$Mismatch_base_in_target_region+=$1;
	}
}
close(TMP);
$Mismatch_rate_in_target_region=$Mismatch_base_in_target_region/$Effective_sequences_on_target;
$Mismatch_rate_in_all_effective_sequence=$Mismatch_base_in_all_effective_sequence/$Total_effective_yield;

$Average_sequencing_depth_on_target=$Effective_sequences_on_target/$Initial_bases_on_target;
$Average_sequencing_depth_near_target=$Effective_sequences_near_target/$Initial_bases_near_target;
	
$Base_covered_on_target=`wc -l $outdir/target.depth | awk '{print \$1}'`;
chomp($Base_covered_on_target);
$Coverage_of_target_region=$Base_covered_on_target/$Initial_bases_on_target;
$Base_covered_near_target=`wc -l $outdir/flanking.depth | awk '{print \$1}'`;
chomp($Base_covered_near_target);
$Coverage_of_flanking_region=$Base_covered_near_target/$Initial_bases_near_target;

my $tmp1=`awk '\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};END{print total1"\t"total2"\t"total3}' $outdir/target.depth`;
chomp($tmp1);
my @info1;
@info1=split /\t/,$tmp1;
if(defined($info1[0]) or $info1[0]=0)
{
$Fraction_of_target_covered_with_at_least_20x=$info1[0]/$Initial_bases_on_target;
}
if(defined($info1[1]) or $info1[1]=0)
{
$Fraction_of_target_covered_with_at_least_10x=$info1[1]/$Initial_bases_on_target;
}
if(defined($info1[2]) or $info1[2]=0)
{
$Fraction_of_target_covered_with_at_least_4x=$info1[2]/$Initial_bases_on_target;
}

my $tmp2=`awk '\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};END{print total1"\t"total2"\t"total3}' $outdir/flanking.depth`;
chomp($tmp2);
my @info2;
@info2=split /\t/,$tmp2;
if(defined($info2[0]) or $info2[0]=0 )
{
$Fraction_of_flanking_region_covered_with_at_least_20x=$info2[0]/$Initial_bases_near_target;
}
if(defined($info2[1]) or $info2[1]=0 )
{
$Fraction_of_flanking_region_covered_with_at_least_10x=$info2[1]/$Initial_bases_near_target;
}
if(defined($info2[2]) or $info2[1]=0 )
{
$Fraction_of_flanking_region_covered_with_at_least_4x=$info2[2]/$Initial_bases_near_target;
}

open STAT,">$outdir/information.xlsx" or die $!;
print STAT "Initial_bases_on_target:\t$Initial_bases_on_target\n";
print STAT "Initial_bases_near_target:\t$Initial_bases_near_target\n";
print STAT "Initial_bases_on_or_near_target:\t$Initial_bases_on_or_near_target\n";
print STAT "Total_effective_reads:\t$Total_effective_reads\n";
printf STAT "Total_effective_yield(Mb):\t%.2f\n",$Total_effective_yield/1000000;
printf STAT "Average_read_length(bp):\t%.2f\n",$Average_read_length;
printf STAT "Effective_sequences_on_target(Mb):\t%.2f\n",$Effective_sequences_on_target/1000000;
printf STAT "Effective_sequences_near_target(Mb):\t%.2f\n",$Effective_sequences_near_target/1000000;
printf STAT "Effective_sequences_on_or_near_target(Mb):\t%.2f\n",$Effective_sequences_on_or_near_target/1000000;
printf STAT "Fraction_of_effective_bases_on_target:\t%.1f%%\n",100*$Fraction_of_effective_bases_on_target;
printf STAT "Fraction_of_effective_bases_on_or_near_target:\t%.1f%%\n",100*$Fraction_of_effective_bases_on_or_near_target;
printf STAT "Average_sequencing_depth_on_target:\t%.2f\n",$Average_sequencing_depth_on_target;
printf STAT "Average_sequencing_depth_near_target:\t%.2f\n",$Average_sequencing_depth_near_target;
printf STAT "Mismatch_rate_in_target_region:\t%.2f%%\n",100*$Mismatch_rate_in_target_region;
printf STAT "Mismatch_rate_in_all_effective_sequence:\t%.2f%%\n",100*$Mismatch_rate_in_all_effective_sequence;
print STAT "Base_covered_on_target:\t$Base_covered_on_target\n";
printf STAT "Coverage_of_target_region:\t%.1f%%\n",100*$Coverage_of_target_region;
print STAT "Base_covered_near_target:\t$Base_covered_near_target\n";
printf STAT "Coverage_of_flanking_region:\t%.1f%%\n",100*$Coverage_of_flanking_region;
printf STAT "Fraction_of_target_covered_with_at_least_20x:\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_20x;
printf STAT "Fraction_of_target_covered_with_at_least_10x:\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_10x;
printf STAT "Fraction_of_target_covered_with_at_least_4x:\t%.1f%%\n",100*$Fraction_of_target_covered_with_at_least_4x;
printf STAT "Fraction_of_flanking_region_covered_with_at_least_20x:\t%.1f%%\n",100*$Fraction_of_flanking_region_covered_with_at_least_20x;
printf STAT "Fraction_of_flanking_region_covered_with_at_least_10x:\t%.1f%%\n",100*$Fraction_of_flanking_region_covered_with_at_least_10x;
printf STAT "Fraction_of_flanking_region_covered_with_at_least_4x:\t%.1f%%\n",100*$Fraction_of_flanking_region_covered_with_at_least_4x;
close STAT;

open DB,"<$outdir/target.depth" or die $!;
open DF,">$outdir/depth_frequency.xls" or die $!;
my %depth=();
while(<DB>)
{
	chomp;
	my @tmp=split;
	$depth{$tmp[2]}++;
}
close(DB);
my $maxCov=0;
$depth{"0"}+=$Initial_bases_on_target-$Base_covered_on_target;
	
foreach my $depth (sort {$a<=>$b} keys %depth)
{
	my $per=$depth{$depth}/$Initial_bases_on_target;
	$maxCov = $per if($per > $maxCov);
	print DF "$depth\t$per\t$depth{$depth}\n";
}
close(DF);
	
open CU,">$outdir/cumu.xls" or die $!;
print CU "Depth\tTRPercent\n";
my @depth= sort {$a<=>$b} keys %depth;

foreach my $depth1 (sort {$a<=>$b} keys %depth)
{
	my $tmp=0;
    foreach my $depth2 (@depth)
    {
		if($depth2 >= $depth1)
        {
        	$tmp+=$depth{$depth2};
        }
    }
    $tmp = $tmp/$Initial_bases_on_target;
    print CU "$depth1\t$tmp\n";
}
close(CU);

if($Plot)
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
	if($Average_sequencing_depth_on_target<30)
	{
		$xlim=100;
		$xbin=20;
	}elsif($Average_sequencing_depth_on_target < 50)
	{
		$xlim=160;
		$xbin=20;
	}elsif($Average_sequencing_depth_on_target < 120)
	{
		$xlim=250;
		$xbin=50;
	}else{
		$xlim=600;
		$xbin=100;
	}
	histPlot($outdir,"$outdir/depth_frequency.xls",$ylim,$ybin,$xlim,$xbin);
	cumuPlot($outdir,"$outdir/cumu.xls",$xlim,$xbin);
}

my %dep=();
my %cov=();
open IN,"$outdir/target.depth" or die "Can't open the $outdir/target.depth:$!";
while (<IN>)
{
	chomp;
	my ($chr,$dd)=(split /\t/,$_)[0,2];
	$cov{$chr}++;
	$dep{$chr}+=$dd;
}
close IN;
open OUT,">$outdir/chrall.stat" or die $!;
print OUT "Chr\tCoverage\tDepth\n";
foreach my $ch (sort keys %hash)
{
	my $cp = 100*$cov{$ch}/$hash{$ch};
	my $dp = $dep{$ch}/$hash{$ch};
	printf OUT "$ch\t%.2f%%\t%.2f\n",$cp,$dp;
}
close OUT;

`gzip -v $outdir/whole_genome.depth $outdir/target.depth $outdir/flanking.depth`;

sub cumuPlot {
        my ($outdir, $dataFile, $xlim, $xbin) = @_;
        my $figFile = "$outdir/cumuPlot.pdf";
        my $Rline=<<Rline;
        pdf(file="$figFile",w=8,h=6)
        rt <- read.table("$dataFile",header=T)
        opar <- par()
        x <- rt\$Depth[1:($xlim+1)]
        y <- 100*rt\$TRPercent[1:($xlim+1)]
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
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
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
        par(opar)
        dev.off()

Rline
        open (ROUT,">$figFile.R");
        print ROUT $Rline;
        close(ROUT);

        system("Rscript $figFile.R");
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
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
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
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);

	system("Rscript $figFile.R");
	#system("rm  $figFile.R  $figFile.Rout");
}
