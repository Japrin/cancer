#!/usr/bin/perl
#============================================================================
# Name        		: ADTEx.zygosity.segments.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sun Dec 28 10:02:01 2014
# Last Modified By	: 
# Last Modified On	: Sun Dec 28 10:02:01 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl ADTEx.zygosity.segments.pl [option] <infile>

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
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

my $preZygosity="";
my $preChr="";
my $preBeg="";
my $preEnd="";
my @preMirroredBAF=();
my $preMarkerN=0;
my $preCN="";
my $preZ="";

open $in,$infile or die "Cann't open file $infile ($!) \n";
print "chr\tbeg\tend\tzygosity\tmean.mirror.baf\tcn\tz\ttarget.number\n";
<$in>;
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    #chrom   SNP_loc control_BAF     tumor_BAF       mirrored_BAF    control_doc     tumor_doc       control_vdoc    tumor_vdoc      cn	z       zygosity
    #1       876499  0.653846        0.222222        0.222222        26      18      17      4       3       1       LOH
    #1       877715  0.45    0.366667        0.366667        20      30      9       11      3       1       LOH
    my ($chr,$pos,$mirroredBAF,$cn,$z,$zygosity)=@F[0,1,4,9,10,11];
    if($chr ne $preChr || $zygosity ne $preZygosity || $cn ne $preCN || $z ne $preZ)
    {
    	# output previous
	if($preChr ne "")
	{
		print "$preChr\t$preBeg\t$preEnd\t$preZygosity\t". (@preMirroredBAF?sum(@preMirroredBAF)/@preMirroredBAF:"NA") ."\t$preCN\t$preZ\t$preMarkerN\n";
	}
	#update 'pre' variables
	$preZygosity=$zygosity;
	$preChr=$chr;
	$preBeg=$pos;
	$preEnd=$pos;
	@preMirroredBAF=($mirroredBAF);
	$preMarkerN=1;
	$preCN=$cn;
	$preZ=$z;
    }else
    {
    	$preEnd=$pos;
	push @preMirroredBAF,$mirroredBAF;
	$preMarkerN++;
    }
}
if($preChr ne "")
{
	print "$preChr\t$preBeg\t$preEnd\t$preZygosity\t". (@preMirroredBAF?sum(@preMirroredBAF)/@preMirroredBAF:"NA") ."\t$preCN\t$preZ\t$preMarkerN\n";
}
###################################################################




sub usage
{
    die `pod2text $0`;
}
