#!/usr/bin/perl
#============================================================================
# Name        		: strandbias_spectrum.pl
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

	perl strandbias_spectrum.pl [option] [infile]

	-o	outDir [default ./]
	-p	plot [default OFF]
	-i	ref column [0-based, default 2] 
	-j	mut column [0-based, default 3]
	-k	strand column [0-based, default 5]
	-s	sampleID [default "this study"]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;

my ($in,$out);
my ($opt_h,$opt_o,$opt_p,$opt_i,$opt_j,$opt_k,$opt_s);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o,"p"=>\$opt_p,"i=i"=>\$opt_i,"j=i"=>\$opt_j,"k=i"=>\$opt_k,"s=s"=>\$opt_s);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_o)) { $opt_o="./"; }
if(!defined($opt_i)) { $opt_i=2; }
if(!defined($opt_j)) { $opt_j=3; }
if(!defined($opt_k)) { $opt_k=5; }
if(!defined($opt_s)) { $opt_s="this study"; }
mkpath $opt_o;

my %count=(
	"T"=>{
		'TG'=>0,
		'TC'=>0,
		'TA'=>0,
		'CA'=>0,
		'CT'=>0,
		'CG'=>0,
	},
	"UT"=>{
		'TG'=>0,
		'TC'=>0,
		'TA'=>0,
		'CA'=>0,
		'CT'=>0,
		'CG'=>0,
	}
);
while(<>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	#chr1    17071876        .       C       T
	my @F=split /\t/;
	my ($from,$to,$strand)=@F[$opt_i,$opt_j,$opt_k];

	my $mutStrand="";
	if($from ne "T" && $from ne "C")
	{
		$from=~tr/ATCG/TAGC/;
		$to=~tr/ATCG/TAGC/;
		$mutStrand="-";
	}else
	{
		$mutStrand="+";
	}
	if($mutStrand eq $strand) { $mutStrand="T"; }
	else{ $mutStrand="UT"; }
	$count{$mutStrand}->{"$from$to"}++;
}

my @title=("T>G","T>C","T>A","C>A","C>T","C>G");
open $out,">","$opt_o/$opt_s.strandbias_spectrum.txt" or die "$!";
printf $out "\t%s\n",join("\t",@title);
printf $out "T";
foreach ('TG', 'TC', 'TA', 'CA', 'CT', 'CG') { printf $out "\t%s",$count{'T'}->{$_}; }
printf $out "\n";
printf $out "UT";
foreach ('TG', 'TC', 'TA', 'CA', 'CT', 'CG') { printf $out "\t%s",$count{'UT'}->{$_}; }
printf $out "\n";
close $out;

if(!defined($opt_p)) { exit 0; }

my $RCommand=sprintf <<Here;
require(plotrix)
mybarplot <- function (x, ...) 
{ 
	ya<-(ceiling(max(x)/5))*7
	at<-barplot(x,col=c("darkred","skyblue"),space=c(0,0.1),beside=T,ylab="Number of mutations",axisnames=F,main="Strand bias of Mutation ($opt_s)",cex.axis=1.5,cex.lab=1.5,cex.main=1.5,ylim=c(0,ya),...) 
	atLab<-apply(at,2,mean)
	staxlab(1,atLab,colnames(x),srt=45,cex=1.3)	
	text(at,x+ya/33,labels=as.character(x))
	legend("top",legend=c("Transcribed","Untranscribed"),bty="n",horiz=T,fill=c("darkred","skyblue"),cex=2)
}
setwd(\"$opt_o\")
a<-read.table("$opt_s.strandbias_spectrum.txt",header=T,sep="",row.names=1,check.names=F)
a.matrix=as.matrix(a)
print(a.matrix)
n<-dim(a.matrix)[1]
row.sum<-apply(a.matrix,1,sum)
#a.matrix<-a.matrix/row.sum
print(row.sum)

png("$opt_s.Strandbias_spectrum.png",width=800,height=600)
par(mar=c(5,5,4,2)+0.1,xpd=F)
mybarplot(a.matrix)
dev.off()
Here
open $out,">","$opt_o/strandbias_spectrum.R" or die "$!";
print $out "$RCommand";
close $out;
system("/usr/local/bin/R CMD BATCH $opt_o/strandbias_spectrum.R >$opt_o/strandbias_spectrum.Rout  2> $opt_o/strandbias_spectrum.Rout");


############################################################################
sub usage
{
	die `pod2text $0`;
}
