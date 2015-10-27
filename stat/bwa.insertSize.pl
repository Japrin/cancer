#!/usr/bin/perl
#============================================================================
# Name        		: bwa.insertSize.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon Jul 22 13:00:22 2013
# Last Modified By	: 
# Last Modified On	: Mon Jul 22 13:00:22 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl bwa.insertSize.pl [option] <inbam> <output prefix>

	-q	mapping quality must be largert than INT [default 30]
	-l	number of alignment used to evaluation; if not specified all; if specified, must be larger than 10000
	-isize	maximum isize in the plot [default 500]
	-sampleID	sampleID [default ""]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use File::Path;

my ($in,$out);
my ($opt_h,$opt_q,$opt_l,$opt_isize,$opt_sampleID);
GetOptions("h"	=>\$opt_h,"q=i"=>\$opt_q,"l=i"=>\$opt_l,"isize=i"=>\$opt_isize,"sampleID=s"=>\$opt_sampleID);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_q)) { $opt_q=30; }
if(!defined($opt_isize)) { $opt_isize=500; }
if(!defined($opt_sampleID)) { $opt_sampleID=""; }
if(defined($opt_l) && $opt_l<10000) { usage(); }

my $inbam=shift @ARGV;
my $prefix=shift @ARGV;
my $outDir=dirname($prefix);
my $bname=basename($prefix);

my %h=();
open $in, "samtools view $inbam -q $opt_q -F 0x404 |" or die "$!";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	## HWI-ST1276:116:D22F3ACXX:1:2303:10166:8402      pPR1    1       13228   39      100M    =       13426   298     CATCAGGCACCAAAGGGATTCTGCCAGCATAGTGCTCCTGGACCAGTGATACACCCGGCACCCTGTCCTGGACACGCTGTTGGCCTGGATCTGAGCCCTG    CEBFCEEFCFGDDEHHHCBDGFGGGDGGECEHAGFGHHGIHCGHEHCGDEBFFFIH@HHEGHGGGCGIFEEDFEF@HFGEDGHHHGGGDDHGGDFGFFDD    X0:i:1  X1:i:1  XA:Z:15,-102517843,100M,1;      MD:Z:100        RG:Z:656.1      XG:i:0  AM:i:16 NM:i:0  SM:i:23 XM:i:0  XO:i:0  MQ:i:39 OQ:Z:CCCFFFFFHHHHHJJJJIIJJJIJJIIJJHJJGIIJJJJJJIJJJJGHJJIIIIJJJIJJHHHFFFFFEEEEEDDDBBDBDDDDDDDDDDDDDCDDDDDD       XT:A:U
	my ($id,$isize)=@F[0,8];
	$h{$id}=abs($isize);
	if(defined($opt_l))
	{
		my $n=keys %h;
		if($n>$opt_l) { last; }
	}
}
#my %isize=();
#foreach (values %h) { $isize{$_}++; }
#open $out,">","$prefix.bwa.insertSize.list" or die "$!";
#for(my $i=0;$i<=$opt_isize;$i++) 
#{
#	printf $out "$i\t%d\n",$isize{$i}?$isize{$i}:0;
#}
#close $out;
open $out,">","$prefix.bwa.insertSize.list" or die "$!";
foreach (values %h) { print $out "$_\n"; }
close $out;
my $RCommand=sprintf <<Here;
library(RColorBrewer)
setwd(\"$outDir\")
a<-read.table("$bname.bwa.insertSize.list",header=F,sep="\\t")
names(a)<-c("isize")
f<-a\$isize<5000 & a\$isize>0
a.hist<-hist(a\$isize[f],breaks=500,plot=F)
print(a.hist)
png("$bname.bwa.insertSize.png",width=800,height=600)
#hist(a\$isize[f],prob=T,breaks=500,xlim=c(0,500),col="red",pch=20,xlab="Insert Size",ylab="Density",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Insert Size Of $opt_sampleID")
#lines(density(a\$isize[f]),xlim=c(0,500),col="darkgreen",lw=4)
plot(density(a\$isize[f]),xlim=c(0,500),col="darkgreen",type="l",lwd=4,xlab="Insert Size",ylab="Density",cex.lab=1.5,cex.axis=1.5,cex.main=1.5,main="Insert Size Of $opt_sampleID")
dev.off()
Here
open $out,">","$prefix.bwa.insertSize.R" or die "$!";
print $out "$RCommand";
close $out;
system("/usr/local/bin/R CMD BATCH $prefix.bwa.insertSize.R");

###################################################################



sub readList
{
	my $in;
	my ($pList,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
	}
}

sub usage
{
	die `pod2text $0`;
}
