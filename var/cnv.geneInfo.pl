#!/usr/bin/perl
#============================================================================
# Name        		: cnv.geneInfo.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Wed Dec 25 11:20:26 2013
# Last Modified By	: 
# Last Modified On	: Wed Dec 25 11:20:26 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl cnv.geneInfo.pl [option] <infile>
	-f	gene field [default "Gene"]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_f);
GetOptions("h"	=>\$opt_h,"f=s"=>\$opt_f);
if(@ARGV<1 || $opt_h) { usage(); }
my $infile=shift @ARGV;
if(!defined($opt_f)) { $opt_f="Gene"; }

open $in,$infile or die "Cann't open file $infile ($!) \n";
$_=<$in>;
chomp $_;
my @F=split /\t/,$_;
my $iGene="";
my $iType="";
my @_setNa=();
for(my $i=0;$i<@F;$i++)
{
	if($F[$i] eq "$opt_f")
	{
		$iGene=$i;
	}
	if($F[$i] eq "SVType")
	{
		$iType=$i;
	}
	if($F[$i] eq "Func" || $F[$i] eq "ExonicFunc" || $F[$i] eq "AAChange" || $F[$i] eq "cpgIslandExt" || $F[$i] eq "evofold" || $F[$i] eq "wgRna" || $F[$i] eq "targetScanS" || $F[$i] eq "phastConsElements46way" || $F[$i] eq "genomicSuperDups" || $F[$i] eq "tfbsConsSites" || $F[$i] eq "cytoband" || $F[$i] eq "Repeat")
	{
		push @_setNa,$i;
	}
}
print "$_\n";
#print STDERR "$_\n";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#Chr     Start   End     Ref     Alt     Func.knownGene  Gene.knownGene  ExonicFunc.knownGene    AAChange.knownGene      cpgIslandExt    evofold wgRna   targetScanS     phastConsElements46way  genomicSuperDups        tfbsConsSites   cytoband        Repeat  CopyNumber      Size    SVID    SVType
	#1       72051936        72305640        0       0       exonic  NEGR1   .       .       .       Score=98;Name=82427_0_+_98      .       .       Score=725;Name=lod=1184 .       Score=1000;Name=V$POU3F2_02,V$AML1_01,V$SRY_01,V$NKX25_02       1p31.1  Score=20688;Name="157703:L1PA10(LINE)"  1       253704  1       loss
	#1       72051936        72051936        0       0       intronic        NEGR1   .       .       .       .       .       .       .       .       .       1p31.1  Score=361;Name="157323:MIRc(SINE)"      1       253704  1       breakpoint
	#1       72305640        72305640        0       0       intronic        NEGR1   .       .       .       .       .       .       Score=381;Name=lod=47   .       .       1p31.1  .       1 #       253704  1       breakpoint
	#if($F[$iType] eq "breakpoint")
	#{
	#	#print "$_\n";
	#}else
	{
		my @g=split /[,;]/,$F[$iGene];
		my @_o=@F;
		#print "$line\n";
		foreach (@g)
		{
			#if(/dist=/) { next; }
			$_o[$iGene]=$_;
			for(my $_j=0;$_j<@_setNa;$_j++)
			{
				$_o[$_setNa[$_j]]=".";
			}
			print join("\t",@_o)."\n";
		}
	}
	
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
