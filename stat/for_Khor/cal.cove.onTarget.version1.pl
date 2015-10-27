#!/usr/bin/perl
=pod

=head1 Usage

	perl cal.cove.onTarget.pl [option] <bam> <bed> <genome>

	-o	output prefix
	-e	flank region [default 100]
	-h	display this help and exit

	bam: input bam file
	bed: target region file, bed format
	genome: 
		<chromName><TAB><chromSize>
		For example, Human (hg19):
		chr1	249250621
		chr2	243199373

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;
use File::Basename;

my ($in,$out);
my ($opt_h,$opt_o,$opt_e);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o,"e=i"=>\$opt_e);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_o)) { $opt_o="."; }
if(!defined($opt_e)) { $opt_e=100; }

my $bamFile=shift @ARGV;
my $bedFile=shift @ARGV;
my $genomeFile=shift @ARGV;

my $bedBase=basename($bedFile);
$bedBase=~s/\.bed$//;

mkpath($opt_o);
`sort -k 1,1 -k 2g,2 -k 3g,3 $bedFile > $opt_o/$bedBase.sort.bed`;
#`slopBed -i $opt_o/$bedBase.sort.bed -b $opt_e -g $genomeFile | mergeBed -i stdin > $opt_o/$bedBase.sort.extend${opt_e}.bed`;

my $len_TR=`awk '{sum+=\$3-\$2}END{print sum}' out/TR.sort.bed`;
chomp $len_TR;
my $sequencedAll=`samtools depth $bamFile | awk '{sum+=\$3}END{print sum}'`;
chomp $baseAll;

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

my $onTargetRate=$sequencedOnTarget/$sequencedAll
my $depthOnTarget=$sequencedOnTarget/$coveredOnTarget;
my $coverageRateOnTarget=$coveredOnTarget/$len_TR;

printf "length of target\t$len_TR\n";
printf "sequenced bases\t$sequencedAll\n";
printf "sequenced bases on target\t$sequencedOnTarget\n";
printf "on target rate\t$onTargetRate\n";
printf "average depth on target\t$$depthOnTarget\n";
printf "coverage rate of target\t$coverageRateOnTarget\n";


###################################################################




sub usage
{
	die `pod2text $0`;
}
