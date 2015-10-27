#!/usr/bin/env perl
#============================================================================
# Name        		: snv_spectrum.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sat Mar 16 20:09:26 2013
# Last Modified By	: 
# Last Modified On	: Sat Mar 16 20:09:26 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl snv_spectrum.pl [option] [infile]

	-o	outDir [default ./]
	-p	plot [default OFF]
	-i	ref column [0-based, default 3] 
	-j	mut column [0-based, default 4] 
	-s	sampleID [default "this study"]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;

my ($in,$out);
my ($opt_h,$opt_o,$opt_p,$opt_i,$opt_j,$opt_s);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o,"p"=>\$opt_p,"i=i"=>\$opt_i,"j=i"=>\$opt_j,"s=s"=>\$opt_s);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_o)) { $opt_o="./"; }
if(!defined($opt_i)) { $opt_i=3; }
if(!defined($opt_j)) { $opt_j=4; }
if(!defined($opt_s)) { $opt_s="this study"; }
mkpath $opt_o;

my %spec=(
		'AC'=>0,
		'TG'=>0,
		'AG'=>1,
		'TC'=>1,
		'AT'=>2,
		'TA'=>2,
		'CA'=>3,
		'GT'=>3,
		'CT'=>4,
		'GA'=>4,
		'CG'=>5,
		'GC'=>5,
		'indel'=>6,
);
#my @cate=('A>C,T>G','A>G,T>C','A>T,T>A','C>A,G>T','C>T,G>A','C>G,G>C');
my @cate=("A/T>C/G","A/T>G/C","A/T>T/A","C/G>A/T","C/G>T/A","C/G>G/C","indel","sampleID");
my @count=(0,0,0,0,0,0,0,$opt_s);
while(<>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	#chr1    17071876        .       C       T
	my @F=split /\t/;
	my ($from,$to)=@F[$opt_i,$opt_j];
	$to=(split /,/,$to)[0];
	#printf "$from\t$to\n";
	my $i;
	if(length($from) != length($to))
	{
		$i=$spec{"indel"};
	}else
	{
		$i=$spec{"$from$to"};
	}
	$count[$i]++;
}
open $out,">","$opt_o/spectrum.txt" or die "$!";
printf $out "%s\n",join("\t",@cate);
printf $out "%s\n",join("\t",@count);
close $out;

if(!defined($opt_p)) { exit 0; }

my $RCommand=sprintf <<Here;
require(plotrix)
library(RColorBrewer)
mybarplot <- function (x, ...) 
{ 
	ya<-(floor(max(x)/5))*6
	at<-barplot(x,col=brewer.pal(7,"Set1"),space=c(0,0.1),beside=T,ylab="Number of mutations",axisnames=F,main="Mutation Spectrum Of $opt_s",cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(0,ya),...) 
	at<-apply(at,2,mean)
	staxlab(1,at,colnames(x),srt=45,cex=1.3)
	text(at,x+ya/30,labels=as.character(x))
}
setwd(\"$opt_o\")
a<-read.table("spectrum.txt",header=T,sep="\\t",check.names=F)
a.matrix=as.matrix(a[,1:7])
a.matrix
png("Spectrum.png",width=800,height=600)
par(mar=c(5,5,4,2)+0.1,xpd=T)
mybarplot(a.matrix)
dev.off()
Here
open $out,">","$opt_o/spectrum.R" or die "$!";
print $out "$RCommand";
close $out;
system("R CMD BATCH $opt_o/spectrum.R");


############################################################################
sub usage
{
	die `pod2text $0`;
}
