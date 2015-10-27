#!/usr/bin/perl
#============================================================================
# Name        		: surroundingSequence.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Thu Sep 29 14:51:22 2011
# Last Modified By	: 
# Last Modified On	: Thu Sep 29 14:51:22 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl surroundingSequence.pl [option] <fasta> <infile>

	-o	outDir [default ./]
	-p	plot [default OFF]
	-t	vcf or bed [default vcf]
	-s	surrounding size [defualt 10]
	-triplet output 'triplet' file [default OFF]
	-sample	sample id [default ""]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Bio::DB::Sam;
use Data::Dumper;
use File::Path;

my ($in,$out);
my ($opt_h,$opt_t,$opt_s,$opt_o,$opt_p,$opt_triplet,$opt_sample);
GetOptions("h"	=>\$opt_h,"t=s"=>\$opt_t,"s=i"=>\$opt_s,"o=s"=>\$opt_o,"p"=>\$opt_p,"triplet"=>\$opt_triplet,"sample=s"=>\$opt_sample);
if(@ARGV<2 || $opt_h) { usage(); }
my $fafile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_t)) { $opt_t="vcf"; }
if($opt_t ne "vcf" && $opt_t ne "bed") { usage(); }
if(!defined($opt_s)) { $opt_s=10; }
if(!defined($opt_o)) { $opt_o="./"; }
if(!defined($opt_sample)) { $opt_sample=""; }
mkpath $opt_o;

#5 prime and 3prime flank bases
my %mutMap=("AA"=>0,"AT"=>1,"AC"=>2,"AG"=>3,"TA"=>4,"TT"=>5,"TC"=>6,"TG"=>7,"CA"=>8,"CT"=>9,"CC"=>10,"CG"=>11,"GA"=>12,"GT"=>13,"GC"=>14,"GG"=>15);
my %Triplet=(
		'CA'=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'CT'=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'CG'=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'TA'=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'TC'=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		'TG'=>[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
	);

if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "$!"; }
else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
my $fai = Bio::DB::Sam::Fai->load($fafile);

my %context=('all'=>[],"A"=>[],"T"=>[],"C"=>[],"G"=>[]);
for(my $i=0;$i<$opt_s*2+1;$i++)
{
	push @{$context{'all'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'TG'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'TC'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'TA'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'CA'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'CT'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
	push @{$context{'CG'}},{"A"=>0,"T"=>0,"C"=>0,"G"=>0};
}
#open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Adenoma
	my @field=split /\t/;
	my ($chr,$pos,$ref,$alt);
	if($opt_t eq "vcf")
	{
		($chr,$pos,$ref,$alt)=@field[0,1,3,4];
	}elsif($opt_t eq "bed")
	{
		($chr,$pos,$ref,$alt)=@field[0,2,3,4];
	}
	
	my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-$opt_s,$pos+$opt_s));
	$dna_string="\U$dna_string";
	
	if($ref ne "T" && $ref ne "C")
	{
		$ref=~tr/ATCG/TAGC/;
		$alt=~tr/ATCG/TAGC/;
		$dna_string=~tr/ATCG/TAGC/;
		$dna_string=reverse($dna_string);
	}
	
	my @ss=split //,$dna_string;
	for(my $i=0;$i<@ss;$i++)
	{
		$context{'all'}->[$i]->{$ss[$i]}++;
		$context{"$ref$alt"}->[$i]->{$ss[$i]}++;
	}
	#printf "$chr\t%d\t%d\t%s\n",$pos-1,$pos,$dna_string;
	if($opt_triplet)
	{
		my $_base5=substr($dna_string,$opt_s-1,1);
		my $_base3=substr($dna_string,$opt_s+1,1);
		$Triplet{"$ref$alt"}->[$mutMap{"$_base5$_base3"}]++;
	}
}
if($opt_triplet)
{
	open $out,">","$opt_o/$opt_sample.mutation_context_triplet.txt" or die "$!";
	my @_c=sort keys %mutMap;
	printf $out "\t%s\n",join("\t",@_c);
	foreach my $_m (keys %Triplet)
	{
		printf $out "$_m";
		for(my $i=0;$i<@_c;$i++)
		{
			printf $out "\t%d",$Triplet{$_m}->[$i];
		}
		printf $out "\n";
	}
	close $out;
}
foreach my $k ("TG","TC","TA","CA","CT","CG")
{
	open $out,">","$opt_o/$opt_sample.mutation_context_$k.txt" or die "$!";
	#printf "mutation context of $k:\n";
	for(my $_ii=-$opt_s;$_ii<=$opt_s;$_ii++) { printf $out "\t%d",$_ii; }
	printf $out "\n";
	my $p=$context{$k};
	foreach my $i ("A","T","C","G")
	{
		printf $out "$i";
		for(my $j=0;$j<@$p;$j++)
		{
			printf $out "\t%s",$p->[$j]->{$i};
		}
		printf $out "\n";
	}
	close $out;
	if(!defined($opt_p)) { next; }

	my ($from,$to)=split //,$k;
	my $title="$from>$to";
	my $RCommand=sprintf <<Here;
library(RColorBrewer)
setwd(\"$opt_o\")
a<-read.table("$opt_sample.mutation_context_$k.txt",header=T,sep="\\t",row.names=1,check.names=F)
a.matrix=as.matrix(a)
n<-sum(a.matrix[,1])
b.matrix<-a.matrix/n
png("$opt_sample.mutation_context_$k.png",width=800,height=600)
par(mar=c(5,5,4,4)+0.1,xpd=T)
barplot(b.matrix,col=brewer.pal(4,"Set1"),border=FALSE,space=c(0),xlab="Postion with respect to mutated base",ylab="Percentage of bases",main=paste("mutation context of $title (n=",n,")",sep=""),cex.axis=1.5,cex.lab=1.5,cex.main=1.5,cex.names=1.5)
legend("right",legend=row.names(b.matrix),bty="n",horiz=F,fill=brewer.pal(4,"Set1"),inset=c(-0.09,0),cex=2)
dev.off()
Here
	open $out,">","$opt_o/$opt_sample.mutation_context_$k.R" or die "$!";
	print $out "$RCommand";
	close $out;
	system("R CMD BATCH $opt_o/$opt_sample.mutation_context_$k.R");
}




############################################################################
sub usage
{
	die `pod2text $0`;
}
