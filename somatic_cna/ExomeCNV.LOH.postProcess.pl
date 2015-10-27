#!/usr/bin/perl
#============================================================================
# Name        		: ExomeCNV.LOH.postProcess.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Fri Dec 26 14:21:19 2014
# Last Modified By	: 
# Last Modified On	: Fri Dec 26 14:21:19 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl ExomeCNV.LOH.postProcess.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use List::Util qw(sum);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }

my $preRegion="";
my $preSegInfo="";
my @mirrorBAF=();
my @mirrorLOHBAF=();
my @mirrorBAFNormal=();
my $marker_number=0;
my $marker_LOH_number=0;

print "chr\tbeg\tend\tnormal_t\tnormal_v\ttumor_t\ttumor_v\tisLOH\tpValue\ttarget.number\ttarget.LOH.number\tmean.mirror.baf\tmean.LOH.mirror.baf\tmean.mirror.baf.normal\n";
while(<>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    #chr5	163260	52235583	31233	14791	43950	22818	TRUE	0.000104283058297967	chr5	163260	163261	179	337	84	461	TRUE	4.92253677832546e-25	1
    #chr5	163260	52235583	31233	14791	43950	22818	TRUE	0.000104283058297967	chr5	413594	413595	25	63	81	102	TRUE	3.48646610840405e-07	1
    my ($chr,$beg,$end,$m_normal_v,$m_normal_t,$m_tumor_v,$m_tumor_t,$m_is_LOH)=@F[0,1,2,12,13,14,15,16];
    my $curRegion="$chr:$beg-$end";
    my $curSegInfo=join("\t",@F[0..8]);
    if($curRegion ne $preRegion)
    {
    	## process pre
	if($preRegion ne "")
	{
		print "$preSegInfo\t$marker_number\t$marker_LOH_number\t". (@mirrorBAF ? sum(@mirrorBAF)/@mirrorBAF : "NA") ."\t". (@mirrorLOHBAF ? sum(@mirrorLOHBAF)/@mirrorLOHBAF : "NA") ."\t". (@mirrorBAFNormal ? sum(@mirrorBAFNormal)/@mirrorBAFNormal : "NA") ."\n";
	}
	## update pre
	$preRegion=$curRegion;
	$preSegInfo=$curSegInfo;
	@mirrorBAF=();
	@mirrorLOHBAF=();
	@mirrorBAFNormal=();
	$marker_number=0;
	$marker_LOH_number=0;
    }
    my $m_baf=$m_tumor_v/$m_tumor_t;
    if($m_baf<0.5) { $m_baf=1-$m_baf; }
    push @mirrorBAF,$m_baf;
    $marker_number++;
    my $m_baf_normal=$m_normal_v/$m_normal_t;
    if($m_baf_normal<0.5) { $m_baf_normal=1-$m_baf_normal; }
    push @mirrorBAFNormal,$m_baf_normal;

    if($m_is_LOH eq "TRUE") 
    { 
    	$marker_LOH_number++;
	push @mirrorLOHBAF,$m_baf;
    }
}
if($preRegion ne "")
{
	print "$preSegInfo\t$marker_number\t$marker_LOH_number\t". (@mirrorBAF ? sum(@mirrorBAF)/@mirrorBAF : "NA") ."\t". (@mirrorLOHBAF ? sum(@mirrorLOHBAF)/@mirrorLOHBAF : "NA") ."\t". (@mirrorBAFNormal ? sum(@mirrorBAFNormal)/@mirrorBAFNormal : "NA") ."\n";
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
