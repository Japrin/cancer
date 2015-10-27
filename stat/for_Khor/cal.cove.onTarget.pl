#!/usr/bin/perl
=pod

=head1 Usage

	perl cal.cove.onTarget.pl [option] <bam> <bed>

	-o	output prefix
	-e	flank region [default 100]
	-h	display this help and exit

	bam: input bam file
	bed: target region file, bed format
	
=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;
use File::Basename;
use Number::Format;

my ($in,$out);
my ($opt_h,$opt_o,$opt_e);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o,"e=i"=>\$opt_e);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_o)) { $opt_o="."; }
if(!defined($opt_e)) { $opt_e=100; }

my $bamFile=shift @ARGV;
my $bedFile=shift @ARGV;

my $bedBase=basename($bedFile);
$bedBase=~s/\.bed$//;

mkpath($opt_o);
`sort -k 1,1 -k 2g,2 -k 3g,3 $bedFile > $opt_o/$bedBase.sort.bed`;

my $len_TR=`awk '{sum+=\$3-\$2}END{print sum}' $opt_o/$bedBase.sort.bed`;
chomp $len_TR;
my $sequencedAll=0;
my $coveredAll=0;
my $t=`samtools depth $bamFile | awk -v OFS="\\t" '{sum+=\$3;cov++;}END{print sum,cov}'`;
chomp $t;
($sequencedAll,$coveredAll)=split /\t/,$t;

open $in,"samtools depth -b $opt_o/$bedBase.sort.bed $bamFile |" or die "$!";
my $sequencedOnTarget=0;
my $coveredOnTarget=0;
while(<$in>)
{
	chomp $_;
	my @F=split /\t/,$_;
	$coveredOnTarget++;
	$sequencedOnTarget+=$F[2];
}

my $onTargetRate=$sequencedOnTarget/$sequencedAll;
my $depthOnTarget=$sequencedOnTarget/$coveredOnTarget;
my $coverageRateOnTarget=$coveredOnTarget/$len_TR;

my $aFormat = new Number::Format;
printf "length of target\t%s\n",$aFormat->format_number($len_TR);
printf "sequenced bases\t%s\n",$aFormat->format_number($sequencedAll);
printf "sequenced bases on target\t%s\n",$aFormat->format_number($sequencedOnTarget);
printf "on target rate\t%4.2f%%\n",$onTargetRate*100;
printf "covered \t%s\n",$aFormat->format_number($coveredAll);
printf "covered on target\t%s\n",$aFormat->format_number($coveredOnTarget);
printf "average depth on target\t%4.2f\n",$depthOnTarget;
printf "coverage rate of target\t%4.2f%%\n",$coverageRateOnTarget*100;


###################################################################




sub usage
{
	die `pod2text $0`;
}
