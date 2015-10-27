#!/usr/bin/perl
#============================================================================
# Name        		: intermutaion.distance.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sat Apr 12 18:46:34 2014
# Last Modified By	: 
# Last Modified On	: Sat Apr 12 18:46:34 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl intermutaion.distance.pl [option] <infile>

	-o	output prefix [required]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;

my ($in,$out);
my ($opt_h,$opt_o);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o);
if(@ARGV<1 || $opt_h) { usage(); }
if(!defined($opt_o)) { usage(); }
my $infile=shift @ARGV;

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
my @cate=("A/T>C/G","A/T>G/C","A/T>T/A","C/G>A/T","C/G>T/A","C/G>G/C","indel");

my %m=();
if($infile=~/\.gz$/){ open $in,"bgzip -cd $infile |" or die "$!"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	my ($chr,$pos)=@F[0,1];
	$chr=~s/^chr//;
	my $ref=$F[3];
	my $alt=(split /,/,$F[4])[0];
	if(!defined($m{$chr}))
	{
		$m{$chr}={};
	}
	$m{$chr}->{$pos}=$spec{"$ref$alt"};
}
my $x=0;
open $out,">","$opt_o.intermutation_distance.txt" or die "$!";
print $out "chr\tpos\tmutType\tx\tdistance\n";
#foreach my $chr (keys %m)
foreach my $chr (1..22,"X")
{
	my @Ary=sort {$a<=>$b} keys %{$m{$chr}};
	my $n=scalar @Ary;
	my $d=0;
	if($n>1) { $d=$Ary[1]-$Ary[0]; }
	$x++;
	printf $out "$chr\t%s\t%s\t%s\t%s\n",$Ary[0],$cate[$m{$chr}->{$Ary[0]}],$x,$d;
	for(my $i=1;$i<$n;$i++)
	{
		$x++;
		printf $out "$chr\t%s\t%s\t%s\t%s\n",$Ary[$i],$cate[$m{$chr}->{$Ary[$i]}],$x,$Ary[$i]-$Ary[$i-1];
	}
}
plot($opt_o);
###################################################################

sub plot
{
	my $output_prefix=shift @_;
	my $infile="$output_prefix.intermutation_distance.txt";
	my $outfile="$output_prefix.intermutation_distance.png";
	my $outR="$output_prefix.intermutation_distance.R";
	my $wDir=dirname($infile);
	my $RCommand=sprintf <<Here;
setwd(\"$wDir\")
library(RColorBrewer)
pal<-brewer.pal(6,"Set1")
##pal<-c("darkblue","black","darkred","purple","yellow","green")
a<-read.table("$infile",header=T)
x_max<-max(a\$x)
y_max<-max(log10(a\$distance))
png("$outfile",width=1600,height=600)
par(mar=c(6,5,2,10.5)+0.1)
plot(x=a\$x,y=a\$distance,type="n",xlim=c(0,x_max),ylim=c(0,y_max),xlab="mutation number",ylab="intermutation distance",main="",yaxt="n",cex.lab=2,cex.axis=1.5)
mType<-c("A/T>C/G","A/T>G/C","A/T>T/A","C/G>A/T","C/G>T/A","C/G>G/C")
for(i in 1:length(mType))
{
	bbb<-a[which(a\$mutType==mType[i]),]
	points(bbb\$x,log10(bbb\$distance),pch=16,cex=1,col=pal[i])
}
axis(2,at=1:as.integer(y_max-0.5+1),labels=10^(1:as.integer(y_max-0.5+1)),cex.axis=1.5)
legend("right",legend=mType,fill=pal,inset=-0.1,horiz=F,xpd=T,cex=1.5)
dev.off()
Here
	open $out,">","$outR" or die "$!";
	print $out "$RCommand";
	close $out;
	system("R CMD BATCH $outR");


}


sub usage
{
	die `pod2text $0`;
}
