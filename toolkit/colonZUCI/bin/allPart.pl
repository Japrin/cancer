#!/usr/bin/perl
#============================================================================
# Name        		: allPart.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Sun Oct 16 20:39:06 2011
# Last Modified By	: 
# Last Modified On	: Sun Oct 16 20:39:06 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl allPart.pl [option] 
	-a	whether output alignment info
	-b	file with another column at the beginning
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_a,$opt_b);
GetOptions("h"	=>\$opt_h,"a"=>\$opt_a,"b"=>\$opt_b);
if($opt_h) { usage(); }

while(<>)
{
	chomp;
	if(/^\s*$/ || /^#/) { next; }
	my @field=split /\t/;
	my $k="";
	if($opt_b) { $k=shift @field; }
	my ($chr,$pos,$rsID,$ref,$alt,$qual,$filter)=@field[0..6];
#ML_ALT=G;TumorRefDepth=10;Tumor_ML_ALT_Dep=4;NormalRefDepth=21;Normal_ML_ALT_Dep=0;NormalRefBQMean=26.71;NormalAltBQMean=0.00;TumorRefBQMean=27.90;TumorAltBQMean=21.25;NormalRefMQMean=51.72;NormalAltMQMean=0.00;TumorRefMQMean=37.00;TumorAltMQMean=34.39;NormalMQ0=0;TumorMQ0=0;NormalIndel=0;TumorIndel=0;DIST_MEAN=4.50;NormalStrandBias=1.0000;TumorStrandBias=1.0000;NormalBQBias=-1.0000;TumorBQBias=0.0000;NormalMQBias=-1.0000;TumorMQBias=0.3924;frequencyP=0.019118;
	my ($Func,$Gene,$ExonicFunc,$AAChange,$Conserved,$SegDup,$SIFT,$PolyPhen2,$LJB_PhyloP,$LJB_MutationTaster,$LJB_LRT,$validation);
	my $info="";
	/Func=(.+?);/ && ($info.="Func=$1;");
	/Gene=(.+?);/ && ($info.="Gene=$1;");
	/ExonicFunc=(.+?);/ && ($info.="ExonicFunc=$1;");
	/AAChange=(.+?);/ && ($info.="AAChange=$1;");
	/Conserved=(.+?);/ && ($info.="Conserved=$1;");
	/SegDup=(.+?);/ && ($info.="SegDup=$1;");
	/tfbs=(.+?);/ && ($info.="tfbs=$1;");
	/miRNA=(.+?);/ && ($info.="miRNA=$1;");
	/SIFT=(.+?);/ && ($info.="SIFT=$1;");
	/PolyPhen2=(.+?);/ && ($info.="PolyPhen2=$1;");
	/LJB_PhyloP=(.+?);/ && ($info.="LJB_PhyloP=$1;");
	/LJB_MutationTaster=(.+?);/ && ($info.="LJB_MutationTaster=$1;");
	/LJB_LRT=(.+?);/ && ($info.="LJB_LRT=$1;");
	/validation=(.+?);/ && ($info.="validation=$1;");
	if($opt_a)
	{
		/ML_ALT=(.+?);/ && ($info.="ML_ALT=$1;");
		/NormalRefDepth=(.+?);/ && ($info.="NormalRefDepth=$1;");
		/Normal_ML_ALT_Dep=(.+?);/ && ($info.="Normal_ML_ALT_Dep=$1;");
		/TumorRefDepth=(.+?);/ && ($info.="TumorRefDepth=$1;");
		/Tumor_ML_ALT_Dep=(.+?);/ && ($info.="Tumor_ML_ALT_Dep=$1;");
		/NormalMQ0=(.+?);/ && ($info.="NormalMQ0=$1;");
		/TumorMQ0=(.+?);/ && ($info.="TumorMQ0=$1;");
		/NormalIndel=(.+?);/ && ($info.="NormalIndel=$1;");
		/TumorIndel=(.+?);/ && ($info.="TumorIndel=$1;");
		/NormalRefBQMean=(.+?);/ && ($info.="NormalRefBQMean=$1;");
		/NormalAltBQMean=(.+?);/ && ($info.="NormalAltBQMean=$1;");
		/TumorRefBQMean=(.+?);/ && ($info.="TumorRefBQMean=$1;");
		/TumorAltBQMean=(.+?);/ && ($info.="TumorAltBQMean=$1;");
		/NormalRefMQMean=(.+?);/ && ($info.="NormalRefMQMean=$1;");
		/NormalAltMQMean=(.+?);/ && ($info.="NormalAltMQMean=$1;");
		/TumorRefMQMean=(.+?);/ && ($info.="TumorRefMQMean=$1;");
		/TumorAltMQMean=(.+?);/ && ($info.="TumorAltMQMean=$1;");
		/NormalStrandBias=(.+?);/ && ($info.="NormalStrandBias=$1;");
		/TumorStrandBias=(.+?);/ && ($info.="TumorStrandBias=$1;");
		/NormalBQBias=(.+?);/ && ($info.="NormalBQBias=$1;");
		/NormalMQBias=(.+?);/ && ($info.="NormalMQBias=$1;");
		/TumorBQBias=(.+?);/ && ($info.="TumorBQBias=$1;");
		/TumorMQBias=(.+?);/ && ($info.="TumorMQBias=$1;");
		/DIST_MEAN=(.+?);/ && ($info.="DIST_MEAN=$1;");
		/frequencyP=(.+?);/ && ($info.="frequencyP=$1;");
	}
	if($k) { $info.="Part=$k;"; }
	printf "%s\t$info\n",join("\t",$chr,$pos,$rsID,$ref,$alt,$qual,$filter);
}



############################################################################
sub usage
{
	die `pod2text $0`;
}
