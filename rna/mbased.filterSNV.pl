#!/usr/bin/perl
#============================================================================
# Name        		: mbased.filterSNV.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Dec  8 17:11:59 2014
# Last Modified By	: 
# Last Modified On	: Mon Dec  8 17:11:59 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl mbased.filterSNV.pl [option] <infile>
    -a	minimum reads for het [default 5]
    -b	minimum reads freq for het [default 0.1]
    -d	minimum snv distance [default 10]
    -n	no functional filter [default OFF]
    -f	filter by "FILTER" [default OFF]
    -g	gender [default "M" ,male;alternative "F",female]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_a,$opt_b,$opt_d,$opt_n,$opt_f,$opt_g);
GetOptions("h"	=>\$opt_h,"a=i"=>\$opt_a,"b=f"=>\$opt_b,"d=i"=>\$opt_d,"n"=>\$opt_n,"f"=>\$opt_f,"g=s"=>\$opt_g);
if(@ARGV<0 || $opt_h) { usage(); }
#my $infile=shift @ARGV;

if(!defined($opt_a)) { $opt_a=5; }
if(!defined($opt_b)) { $opt_b=0.1; }
if(!defined($opt_d)) { $opt_d=10; }
if(!defined($opt_g)) { $opt_g="M"; }

my $hetCalled=0;
my $filteredByReadCount=0;
my $filteredBySNVDist=0;
my $filteredByGene=0;
my $filteredByExon=0;
my $filterByChr=0;
my $nPass=0;

#if($infile=~/\.gz$/){ open $in,"bgzip -cd $infile |" or die "$!"; }
#else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
#while(<$in>)
while(<>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  466.466N
    #GATK:
    #1       877715  rs6605066       C       G       266.77  PASS    Func=intronic;Gene=SAMD11;cpgIslandExt=(Name=CpG:_246);cytoband=1p36.33;1000g2012apr_all=0.89;snp138=rs6605066;ABHom=1.00;AC=2;AF=1.00;AN=2;DP=13;Dels=0.00;FS=0.000;HaplotypeScore=0.5784;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=20.52;VariantType=SNP   GT:AD:DP:GQ:PL  1/1:0,13:13:30:295,30,0
    #1       877782  rs79037098      C       G       404.77  PASS    Func=intronic;Gene=SAMD11;cpgIslandExt=(Name=CpG:_246);cytoband=1p36.33;1000g2012apr_all=0.07;snp138=rs79037098;ABHet=0.333;AC=1;AF=0.500;AN=2;BaseQRankSum=-1.694;DP=30;Dels=0.00;FS=0.000;HaplotypeScore=0.9789;MLEAC=1;MLEAF=0.500;MQ=60.00;MQ0=0;MQRankSum=0.242;QD=13.49;ReadPosRankSum=0.022;VariantType=SNP      GT:AD:DP:GQ:PL  0/1:10,20:30:99:433,0,231
    #
    #SAMTools:
    #chr1    10109   .       A       T       40      .       Func=intergenic;Gene=NONE(dist=NONE),DDX11L1(dist=1765);SegDup=0.99;DP=117;VDB=3.924378e-03;RPB=1.590163e-01;AF1=0.5;AC1=1;DP4=28,35,19,21;MQ=21;FQ=43;PV4=0.84,0.017,0.0031,1  GT:PL:DP:SP:GQ  0/1:70,0,207:103:1:73
    #
    if($opt_f && ($F[6] ne "PASS" && $F[6] ne ".")) { next; }
    
    my ($colFormat,$colSample)=@F[8,9];
    my @colFormat=split /:/,$colFormat;
    my @colSample=split /:/,$colSample;
    my $iGT;
    my $iAD;
    for(my $i=0;$i<@colFormat;$i++)
    {
		if($colFormat[$i] eq "GT") { $iGT=$i; }
		if($colFormat[$i] eq "AD") { $iAD=$i; }
		if(defined($iGT) && defined($iAD)) { last; }
    }
    #if($colSample[$iGT] eq "1/1" || $colSample[$iGT] eq "0/0")
    if($colSample[$iGT] ne "0/1")
    {
    	next;
    }
    $hetCalled++;
   
    my @ad=();
    if($F[7]=~/DP4=(.+?);/)
    {
		#SAMTools
		my $_dp4=$1;
		my @_dp4=split /,/,$_dp4;
		@ad=($_dp4[0]+$_dp4[1],$_dp4[2]+$_dp4[3]);
	
    }else
    {
		#GATK
		@ad=split /,/, $colSample[$iAD];
    }
    if($ad[0] < $opt_a || $ad[1] < $opt_a || $ad[0]/($ad[0]+$ad[1]) < $opt_b || $ad[1]/($ad[0]+$ad[1]) < $opt_b)
    {
        $filteredByReadCount++;
        next;
    }

    my ($chr,$pos,$ref,$alt,$info)=@F[0,1,3,4,7];
    if($info =~ /;SNVDist=(\d+)/ && $1 < $opt_d)
    {
    	$filteredBySNVDist++;
    	next;
    }

    my $gene;
    my @g=();
    if($info =~ /Gene=(.+?);/)
    {
        @g=split /,/,$1;
        $gene=$g[0];
    }

    if(!defined($opt_n))
    {
		if(defined($gene))
		{
			if(@g>1)
			{
				$filteredByGene++;
				next;
			}
		}else
		{
			$filteredByGene++;
			next;
		}

		if($info =~ /Func=(.+?);/)
		{
			if(!($1 eq "exonic" || $1 eq "UTR5" || $1 eq "UTR3"))
			{
				$filteredByExon++;
				next;
			}
		}
    }
    if($chr eq "chrY" || $chr eq "Y")
    {
    	$filterByChr++;
		next;
    }
    if($opt_g eq "M" && ($chr eq "chrX" || $chr eq "X"))
    {
    	$filterByChr++;
		next;
    }
    $nPass++;
    #bed format
    print "$chr\t".($pos-1)."\t$pos\t$ref\t$alt\t$gene\n";
    #print "$line\n";
}

print STDERR "HetCalled\tfilteredByReadCount\tfilteredBySNVDist\tfilteredByGene\tfilteredByExon\tfilterByChromosome\tpassFilter\n";
print STDERR "$hetCalled\t$filteredByReadCount\t$filteredBySNVDist\t$filteredByGene\t$filteredByExon\t$filterByChr\t$nPass\n";

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
