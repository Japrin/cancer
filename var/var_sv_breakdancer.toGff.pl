#!/usr/bin/env perl
#============================================================================
# Name        		: var_sv_breakdancer.toGff.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Wed Dec  4 11:44:15 2013
# Last Modified By	: 
# Last Modified On	: Wed Dec  4 11:44:15 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl var_sv_breakdancer.toGff.pl [option] <infile>

	-id		id [default ""]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_id);
GetOptions("h"	=>\$opt_h,"id=s"=>\$opt_id);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_id)) { $opt_id=""; }

open $in,$infile or die "Cann't open file $infile ($!) \n";
my $i=0;
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	
	#Chr1   Pos1    Orientation1    Chr2    Pos2    Orientation2    Type    Size    Score   num_Reads       num_Reads_lib   Allele_frequency        CH21.final.bam  CH22.final.bam  CH23.final.bam
	#1       10016   16+12-  1       10381   16+12-  ITX     -235    99      2       /PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH22.final.bam|2 -nan    NA      NA      NA
	#1       965794  11+14-  1       966028  1+3-    DEL     151     99      3       /PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH21.final.bam|1:/PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH22.final.bam|1:/PROJ/GR/HUMAN/shidao.cancer.wgs/zhengliangtao/all/aln/CH23.final.bam|1 0.89    0.25    0.19    0.25
	$i++;	
	
	my $size=$F[7];
	my $num_reads=$F[9];
	my $feat=$F[6];

	my $nn=$F[10];
	my @nn=split /:/,$nn;
	for(my $_iii=0;$_iii<@nn;$_iii++)
	{
		my ($_s,$_n)=split /\|/,$nn[$_iii];
		if($_s=~/.+\/(.+?)\./)
		{
			$_s=$1;
		}
		$nn[$_iii]="$_s,$_n";
	}
	$nn=join(",",@nn);

	my $start=($F[1]<=$F[4]?$F[1]:$F[4]);
	my $end=($F[4]>=$F[1]?$F[4]:$F[1]);
	my $tx="";
	
	if ($feat eq "DEL") { $feat="Deletion"; }
	elsif ($feat eq "INS") { $feat="Insertion"; }
	elsif ($feat eq "INV") { $feat="Inversion"; }
	elsif ($feat eq "ITX" || $feat eq "CTX")
	{
		$feat="Translocation";
	}
	my $sv_id="";
	if($opt_id) { $sv_id="$opt_id.$i"; } 
	else { $sv_id="$i"; }

	if($feat eq "Translocation")
	{
		print "$F[0]\tBreakDancer\t$feat\t$F[1]\t$F[1]\t.\t.\t.\tOrientation1=$F[2];Orientation2=$F[5];Score=$F[8];Size=$size;Support=$num_reads;SupportPerID=$nn;TX=$F[6];TCHR=$F[3];TSTART=$F[4];SVID=$sv_id;SVType=$feat\n";
		print "$F[3]\tBreakDancer\t$feat\t$F[4]\t$F[4]\t.\t.\t.\tOrientation1=$F[2];Orientation2=$F[5];Score=$F[8];Size=$size;Support=$num_reads;SupportPerID=$nn;TX=$F[6];TCHR=$F[0];TSTART=$F[1];SVID=$sv_id;SVType=$feat\n";
	}else
	{
		print "$F[0]\tBreakDancer\t$feat\t$start\t$end\t.\t.\t.\tOrientation1=$F[2];Orientation2=$F[5];Score=$F[8];Size=$size;Support=$num_reads;SupportPerID=$nn;TX=na;TCHR=na;TSTART=na;SVID=$sv_id;SVType=$feat\n";
		print "$F[0]\tBreakDancer\t$feat\t$F[1]\t$F[1]\t.\t.\t.\tOrientation1=$F[2];Orientation2=$F[5];Score=$F[8];Size=$size;Support=$num_reads;SupportPerID=$nn;TX=na;TCHR=na;TSTART=na;SVID=$sv_id;SVType=breakpoint\n";
		print "$F[0]\tBreakDancer\t$feat\t$F[4]\t$F[4]\t.\t.\t.\tOrientation1=$F[2];Orientation2=$F[5];Score=$F[8];Size=$size;Support=$num_reads;SupportPerID=$nn;TX=na;TCHR=na;TSTART=na;SVID=$sv_id;SVType=breakpoint\n";
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
