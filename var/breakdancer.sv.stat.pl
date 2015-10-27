#!/usr/bin/env perl
#============================================================================
# Name        		: sv.stat.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Tue Jan 22 16:38:06 2013
# Last Modified By	: 
# Last Modified On	: Tue Jan 22 16:38:06 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl sv.stat.pl [option] <infile>

	-s	sample [default "SAMPLE" ]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_s)) { $opt_s="SAMPLE"; }

my ($iGene,$iType,$iID,$iFunc)=getI($infile);

my %dup=();
my %stat=();
my %weight=(cds=>0,splicing=>1,utr5=>2,utr3=>3,intron=>4,upstream=>5,downstream=>6,ncRNA=>7,intergenic=>8,unknown=>9,length=>10);

open $in,$infile or die "Cann't open file $infile ($!) \n";
$_=<$in>;
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^#/ || /^\s*$/) { next; }
	##chr    start   end     gene_symbol     region  rmsk    segmental duplication   other_info
	#1       1142719 1143137 TNFRSF18        upstream        .       .       Deletion        .       Size 449;Support 45;SVID 1
	#1       1232594 1233098 ACAP3   intronic        .       .       Deletion        .       Size 509;Support 33;SVID 2
	#Chr     Start   End     Ref     Alt     Func    Gene    ExonicFunc      AAChange        evofold wgRna   targetScanS     phastConsElements46way  genomicSuperDups        tfbsConsSites   cytoband        Repeat  Size    Support SupportPerID    TX      TCHR    TSTART  SVID    SVType
	#chr1    1478903 1479004 0       0       UTR3    SSU72   .       .       .       .       .       Score=313;Name=lod=25   .       .       1p36.33 .       -256    6       NHD0815,2,NHD0816,4
	#     na      na      na      1       Insertion
	my @F=split /\t/;
	my $type=$F[$iType];
	if($type eq "breakpoint") { next; }
	if(!exists($stat{$type})) { $stat{$type}={total=>0,cds=>0,splicing=>0,utr5=>0,utr3=>0,intron=>0,upstream=>0,downstream=>0,ncRNA=>0,intergenic=>0,unknown=>0,length=>0}; }
	my $ID=$F[$iID];
	if(defined($ID))
	{
		if($dup{$ID}) { next; }
		$dup{$ID}++;
	}
	my $func=$F[$iFunc];
	my $w=20;
	$stat{$type}->{'total'}++;
	if($func=~/ncRNA/) { $stat{$type}->{'ncRNA'}++; $w=$weight{'ncRNA'};  }
	elsif($func=~/exonic/) { $stat{$type}->{'cds'}++; $w=$weight{'cds'}; }
	elsif($func=~/splicing/) { $stat{$type}->{'splicing'}++; $w=$weight{'splicing'}; }
	elsif($func=~/UTR5/) { $stat{$type}->{'utr5'}++; $w=$weight{'utr5'}; }
	elsif($func=~/UTR3/) { $stat{$type}->{'utr3'}++; $w=$weight{'utr3'}; }
	elsif($func=~/intron/) { $stat{$type}->{'intron'}++; $w=$weight{'intron'}; }
	elsif($func=~/upstream/) { $stat{$type}->{'upstream'}++; $w=$weight{'upstream'}; }
	elsif($func=~/downstream/) { $stat{$type}->{'downstream'}++; $w=$weight{'downstream'}; }
	elsif($func=~/intergenic/) { $stat{$type}->{'intergenic'}++; $w=$weight{'intergenic'}; }
	else { $stat{$type}->{'unknown'}++; $w=$weight{'unknown'}; } 
	if($type!~/Translocation|Insertion|CTX|ITX/) { $stat{$type}->{'length'}+=$F[2]-$F[1]+1; }

}
printf "#sample\tvarType\ttotal\tcds\tsplicing\tutr5\tutr3\tintron\tupstream\tdownstream\tncRNA\tintergenic\tunknown\tlength\n";
foreach my $type (keys %stat)
{
	printf "$opt_s\t$type\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$stat{$type}->{'total'},$stat{$type}->{'cds'},$stat{$type}->{'splicing'},$stat{$type}->{'utr5'},$stat{$type}->{'utr3'},$stat{$type}->{'intron'},$stat{$type}->{'upstream'},$stat{$type}->{'downstream'},$stat{$type}->{'ncRNA'},$stat{$type}->{'intergenic'},$stat{$type}->{'unknown'},$stat{$type}->{'length'};
			
}

############################################################################
sub getI
{
	my ($infile)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	$_=<$in>;
	chomp $_;
	my @F=split /\t/,$_;
	my $iGene="";
	my $iType="";
	my $iID="";
	my $iFunc="";
	for(my $i=0;$i<@F;$i++)
	{
		if($F[$i] eq "Gene")
		{
			$iGene=$i;
		}
		if($F[$i] eq "SVType")
		{
			$iType=$i;
		}
		if($F[$i] eq "SVID")
		{
			$iID=$i;
		}
		if($F[$i] eq "Func")
		{
			$iFunc=$i;
		}
	}
	close $in;
	return ($iGene,$iType,$iID,$iFunc);
}
sub usage
{
	die `pod2text $0`;
}
