#!/usr/bin/perl
#============================================================================
# Name        		: somatic_cnv_varScan.stat.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Fri Oct 11 18:17:16 2013
# Last Modified By	: 
# Last Modified On	: Fri Oct 11 18:17:16 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl somatic_cnv_varScan.stat.pl [option] [infile]

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }

my %stat=();
my $total_count=0;
my $total_size=0;
my $_h=<>;
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#Chr     Start   End     Ref     Alt     Func    Gene    ExonicFunc      AAChange        evofold wgRna   targetScanS     phastConsElements46way  genomicSuperDups        tfbsConsSites   cytoband        Repeat  seg_mean        num_segments    num_markers     event_size size_class      chrom_arm       arm_fraction    chrom_fraction  SVID    SVType
	#1       6592640 6631128 0       0       exonic  NOL9,TAS1R1     .       .       .       .       .       Score=599;Name=lod=362  .  Score=1000;Name=V$AML1_01       1p36.31 Score=11810;Name="9996:HERVFH21-int(LTR)","9992:HERVFH21-int(LTR)","9998:HERVFH21-int(LTR)","9994:HERVFH21-int(LTR)"    -0.744  1       23      38489   focal   1p      0.03%   0.02%   1       deletion 
	#my ($size,$size_class)=$F[-1]=~/event_size=(.+?);size_class=(.+?);/;
	#my $type=$F[7];
	my ($size,$size_class,$type)=@F[20,21,26];
	$stat{$type}->{$size_class}->{'count'}++;
	$stat{$type}->{$size_class}->{'size'}+=$size;
	$total_count++;
	$total_size+=$size;
}
print "CNV Type\tSize Type\tCount\tSize\n";
for my $type ("deletion","amplification")
{
	foreach my $size_class ("focal","large-scale")
	{
		printf "$type\t$size_class\t%d\t%d\n",$stat{$type}->{$size_class}->{'count'}?$stat{$type}->{$size_class}->{'count'}:0,$stat{$type}->{$size_class}->{'size'}?$stat{$type}->{$size_class}->{'size'}:0;
	}
}
printf "Total\t-\t%d\t%d\n",$total_count,$total_size;
###################################################################




sub usage
{
	die `pod2text $0`;
}
