#!/usr/bin/perl
#============================================================================
# Name        		: var.sv.putative.fusion.gene.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Sun Dec 29 17:55:13 2013
# Last Modified By	: 
# Last Modified On	: Sun Dec 29 17:55:13 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl var.sv.putative.fusion.gene.pl [option] <infile>

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

my ($iSVID,$iSVType,$iFunc,$iGene)=getI($infile);

open $in,$infile or die "Cann't open file $infile ($!) \n";
$_=<$in>;
chomp;
print "$_\tFusionType\n";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	my ($partnerOneFunc,$partnerOneGene,$partnerTwoFunc,$partnerTwoGene);
	if($F[$iSVType] eq "DEL" || $F[$iSVType] eq "loss" || $F[$iSVType] eq "Deletion")
	{
		my $lineBP1=<$in>;
		chomp $lineBP1;
		my @_BP1=split /\t/,$lineBP1;
		my $lineBP2=<$in>;
		chomp $lineBP2;
		my @_BP2=split /\t/,$lineBP2;
		if($_BP1[$iFunc] !~ /intergenic|ncRNA/ && $_BP2[$iFunc] !~ /intergenic|ncRNA/)
		{
			if($_BP1[$iGene] ne $_BP2[$iGene])
			{
				print "$line\tDEL_fusion_2_genes\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}elsif($F[$iFunc] eq "exonic" && $_BP1[$iGene] ne "exonic" && $_BP2[$iGene] ne "exonic")
			{
				print "$line\tDEL_exon_inframe\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}
		}
	}elsif($F[$iSVType] eq "INV" || $F[$iSVType] eq "Inversion")
	{
		my $lineBP1=<$in>;
		chomp $lineBP1;
		my @_BP1=split /\t/,$lineBP1;
		my $lineBP2=<$in>;
		chomp $lineBP2;
		my @_BP2=split /\t/,$lineBP2;
		if($_BP1[$iFunc] !~ /intergenic|ncRNA/ && $_BP2[$iFunc] !~ /intergenic|ncRNA/)
		{
			if($_BP1[$iGene] ne $_BP2[$iGene])
			{
				## to do: strand criteria
				print "$line\tINV_fusion_2_genes\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}elsif($F[$iFunc] eq "exonic" && $_BP1[$iGene] ne "exonic" && $_BP2[$iGene] ne "exonic")
			{
				print "$line\tINV_exon_inframe\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}elsif($F[$iFunc] eq "exonic" && ($_BP1[$iGene] eq "exonic" || $_BP2[$iGene] eq "exonic"))
			{
				print "$line\tINV_exon_disrupt\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}

		}

	}elsif($F[$iSVType] eq "INS")  ## tandem duplication
	{
		my $lineBP1=<$in>;
		chomp $lineBP1;
		my @_BP1=split /\t/,$lineBP1;
		my $lineBP2=<$in>;
		chomp $lineBP2;
		my @_BP2=split /\t/,$lineBP2;
		if($_BP1[$iFunc] !~ /intergenic|ncRNA/ && $_BP2[$iFunc] !~ /intergenic|ncRNA/)
		{
			if($_BP1[$iGene] ne $_BP2[$iGene])
			{
				## to do: strand criteria
				print "$line\tTandemDuplication_fusion_2_genes\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}elsif($F[$iFunc] eq "exonic" && $_BP1[$iGene] ne "exonic" && $_BP2[$iGene] ne "exonic")
			{
				print "$line\tINS_exon_inframe\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}elsif($F[$iFunc] eq "exonic" && ($_BP1[$iGene] eq "exonic" || $_BP2[$iGene] eq "exonic"))
			{
				print "$line\tINS_exon_partial\n";
				print "$lineBP1\t.\n";
				print "$lineBP2\t.\n";
			}
			
		}

	}elsif($F[$iSVType] eq "ITX" || $F[$iSVType] eq "CTX" || $F[$iSVType] eq "Translocation")
	{
		my $lineBP1=$line;
		my @_BP1=@F;
		my $lineBP2=<$in>;
		chomp $lineBP2;
		my @_BP2=split /\t/,$lineBP2;
		if($_BP1[$iFunc] !~ /intergenic|ncRNA/ && $_BP2[$iFunc] !~ /intergenic|ncRNA/)
		{
			if($_BP1[$iGene] ne $_BP2[$iGene])
			{
				## to do: strand criteria
				print "$lineBP1\tTranslocation_fusion_2_genes\n";
				print "$lineBP2\tTranslocation_fusion_2_genes\n";
			}else
			{
				print "$lineBP1\tTranslocation_intra_gene\n";
				print "$lineBP2\tTranslocation_intra_gene\n";
			}
		}
	}
	
}
###################################################################

sub getI
{
	my ($infile)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	$_=<$in>;
	chomp;
	my @F=split /\t/;
	my ($iSVID,$iSVType,$iFunc,$iGene);
	for(my $i=0;$i<@F;$i++)
	{
		if($F[$i] eq "SVID")
		{
			$iSVID=$i;
		}elsif($F[$i] eq "SVType")
		{
			$iSVType=$i;
		}elsif($F[$i] eq "Func")
		{
			$iFunc=$i;
		}elsif($F[$i] eq "Gene")
		{
			$iGene=$i;
		}
	}
	return ($iSVID,$iSVType,$iFunc,$iGene);
	close $in;
}


sub usage
{
	die `pod2text $0`;
}
