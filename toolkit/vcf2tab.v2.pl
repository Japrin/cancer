#!/usr/bin/env perl
#============================================================================
# Name        		: vcf2tab.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Aug 15 20:30:46 2011
# Last Modified By	: 
# Last Modified On	: Mon Aug 15 20:30:46 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf2tab.pl [option] 
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if($opt_h) { usage(); }

my @header=();
my %FlagTag=();
my @format=();
my $flag_write=1;
while(<>)
{
	chomp;
	if(/#INFO=<ID=(.+?),.*Type=(.+?),/ &&  $1!~/^(ExcessHet)$/)
	{
		push @header,$1;
		if($2 eq "Flag") { $FlagTag{$1}=1; }
		next;
		#INFO=<ID=Gene,Number=-1,Type=String,Description="Gene name">
	}elsif(/#CHROM/)
	{ 
		@format=split /\t/; 
		splice @format,0,8; 
		next;
	}elsif(/^\s*$/ || /^#/) { next; }
	if($flag_write)
	{
		printf "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\t%s\n",join("\t",@header,@format);
		$flag_write=0;
	}
	##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
	#chr1    1133110 .       G       A       .       .       Func=exonic;Gene=TTLL10;ExonicFunc=synonymous SNV;AAChange=NM_001130045:c.G1905A:p.G635G;ensGene=ENSG00000162571:exonic:synonymous SNV:ENST00000379289:c.G1905A:p.G635G;somaticScore=25;tumorConsensusQuality=33;tumorSNVQuality=33;tumorRMSMappingQuality=60;DepthTumor=2;DepthNormal=1;ML_ALT=A;RefDepth=0;ML_ALT_Dep=2;BQ_MEAN=30.00;MQ_MEAN=60.00;DIST_MEAN=48.50;
	my @field=split /\t/;
	my @out_a=();
	@out_a=@field[0..6];
	my %_t= "$field[7];"=~/(.+?)=(.+?);/g;
	foreach (@header) 
	{ 
		if(!exists($_t{$_})) 
		{ 
			if(exists($FlagTag{$_}))
			{
				if($field[7]=~/$_/){ push @out_a,"Yes"; }
				else{ push @out_a,"No"; }
			}
			else{ push @out_a,"-"; }
			next;
		}
		push @out_a,$_t{$_};
	}
	splice @field,0,8;
	push @out_a,@field;

	printf "%s\n",join("\t",@out_a);

}



############################################################################
sub usage
{
	die `pod2text $0`;
}
