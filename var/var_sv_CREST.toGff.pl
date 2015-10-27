#!/usr/bin/perl
#============================================================================
# Name        		: var_sv_CREST.toGff.pl
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

	perl var_sv_CREST.toGff.pl [option] <infile>

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
	
	#left chr	 left pos	left strand	 # of left soft clipped reads	 right chr	 right pos	 right strand	# right soft clipped reads	 SV type	 coverage at left pos	 coverage at right pos	 assembled length at left_pos	 assembled length at right_pos	 average percent identity at left_pos	 percent of non-unique mapping reads at left_pos	 average percent identity at right_pos	 percent of non-unique mapping reads at right_pos	 start position of consensus mapping to genome	 starting chromosome of consensus mapping	 position of the genomic mapping of consensus starting position	 end position of consensus mapping to genome	 ending chromsome of consnesus mapping	 position of genomic mapping of consensus ending posiiton	consensus sequences
	#chrX    124708914       -       0       chr12   71339066        +       5       CTX     89      17      0       46      0.760867924528302       0.657142857142857       0.98989010989011 0       1       chrX    124708959       140     chr12   71339159        CCCTCAAGTAGGTCCTGGTGTCTGTTCCCTTTTTTGTGTCAATGTGTTCTCATTGTTCAATTCCCACCTACGAGTGAGAACATGCGGTATTTGGTTTTCTGTCCTTGCGATAGTTTGCTGAGAATGATGGTTTCCAGCTT
	$i++;	
	
	
	my $sv_id="";
	if($opt_id) { $sv_id="$opt_id.$i"; } 
	else { $sv_id="$i"; }
	
	my $feat=$F[8];
	my ($annLeft,$annRight)=@F[1,5];
	if($F[0]  eq $F[4] && $annLeft>$annRight)
	{
		($annLeft,$annRight)=($annRight,$annLeft);
	}
	if($feat eq "CTX" || $feat eq "ITX")
	{
		print "$F[0]\tCREST\t$feat\t$annLeft\t$annLeft\t.\t.\t.\tLeftPos=$F[0]:$F[1];RightPos=$F[4]:$F[5];LeftStrand=$F[2];RightStrand=$F[6];LeftSClip=$F[3];RightSClip=$F[7];LeftCoverage=$F[9];RightCoverage=$F[10];LeftAssemblyLength=$F[11];RightAssemblyLength=$F[12];AveragePercentIdentityAtLeftPos=$F[13];PercentOfNonUniqueMappingReadsAtLeftPos=$F[14];AveragePercentIdentityAtRightPos=$F[15];PercentOfNonUniqueMappingReadsAtRightPos=$F[16];StartPositionOfConsensusMappingToGenome=$F[17];StartingChromosomeOfConsensusMapping=$F[18];PositionOfTheGenomicMappingOfConsensusStartingPosition=$F[19];EndPositionOfConsensusMappingToGenome=$F[20];EndingChromsomeOfConsnesusMapping=$F[21];PositionOfGenomicMappingOfConsensusEndingPosiiton=$F[22];ConsensusSequences=$F[23];TCHR=$F[4];TSTART=$F[5];SVID=$sv_id;SVType=$feat\n";
		print "$F[4]\tCREST\t$feat\t$annRight\t$annRight\t.\t.\t.\tLeftPos=$F[0]:$F[1];RightPos=$F[4]:$F[5];LeftStrand=$F[2];RightStrand=$F[6];LeftSClip=$F[3];RightSClip=$F[7];LeftCoverage=$F[9];RightCoverage=$F[10];LeftAssemblyLength=$F[11];RightAssemblyLength=$F[12];AveragePercentIdentityAtLeftPos=$F[13];PercentOfNonUniqueMappingReadsAtLeftPos=$F[14];AveragePercentIdentityAtRightPos=$F[15];PercentOfNonUniqueMappingReadsAtRightPos=$F[16];StartPositionOfConsensusMappingToGenome=$F[17];StartingChromosomeOfConsensusMapping=$F[18];PositionOfTheGenomicMappingOfConsensusStartingPosition=$F[19];EndPositionOfConsensusMappingToGenome=$F[20];EndingChromsomeOfConsnesusMapping=$F[21];PositionOfGenomicMappingOfConsensusEndingPosiiton=$F[22];ConsensusSequences=$F[23];TCHR=$F[0];TSTART=$F[1];SVID=$sv_id;SVType=$feat\n";

	}else
	{

		if($feat eq "INV")
		{
		}
		print "$F[0]\tCREST\t$feat\t$annLeft\t$annRight\t.\t.\t.\tLeftPos=$F[0]:$F[1];RightPos=$F[4]:$F[5];LeftStrand=$F[2];RightStrand=$F[6];LeftSClip=$F[3];RightSClip=$F[7];LeftCoverage=$F[9];RightCoverage=$F[10];LeftAssemblyLength=$F[11];RightAssemblyLength=$F[12];AveragePercentIdentityAtLeftPos=$F[13];PercentOfNonUniqueMappingReadsAtLeftPos=$F[14];AveragePercentIdentityAtRightPos=$F[15];PercentOfNonUniqueMappingReadsAtRightPos=$F[16];StartPositionOfConsensusMappingToGenome=$F[17];StartingChromosomeOfConsensusMapping=$F[18];PositionOfTheGenomicMappingOfConsensusStartingPosition=$F[19];EndPositionOfConsensusMappingToGenome=$F[20];EndingChromsomeOfConsnesusMapping=$F[21];PositionOfGenomicMappingOfConsensusEndingPosiiton=$F[22];ConsensusSequences=$F[23];TCHR=na;TSTART=na;SVID=$sv_id;SVType=$feat\n";
		print "$F[0]\tCREST\t$feat\t$annLeft\t$annLeft\t.\t.\t.\tLeftPos=$F[0]:$F[1];RightPos=$F[4]:$F[5];LeftStrand=$F[2];RightStrand=$F[6];LeftSClip=$F[3];RightSClip=$F[7];LeftCoverage=$F[9];RightCoverage=$F[10];LeftAssemblyLength=$F[11];RightAssemblyLength=$F[12];AveragePercentIdentityAtLeftPos=$F[13];PercentOfNonUniqueMappingReadsAtLeftPos=$F[14];AveragePercentIdentityAtRightPos=$F[15];PercentOfNonUniqueMappingReadsAtRightPos=$F[16];StartPositionOfConsensusMappingToGenome=$F[17];StartingChromosomeOfConsensusMapping=$F[18];PositionOfTheGenomicMappingOfConsensusStartingPosition=$F[19];EndPositionOfConsensusMappingToGenome=$F[20];EndingChromsomeOfConsnesusMapping=$F[21];PositionOfGenomicMappingOfConsensusEndingPosiiton=$F[22];ConsensusSequences=$F[23];TCHR=na;TSTART=na;SVID=$sv_id;SVType=breakpoint\n";
		print "$F[4]\tCREST\t$feat\t$annRight\t$annRight\t.\t.\t.\tLeftPos=$F[0]:$F[1];RightPos=$F[4]:$F[5];LeftStrand=$F[2];RightStrand=$F[6];LeftSClip=$F[3];RightSClip=$F[7];LeftCoverage=$F[9];RightCoverage=$F[10];LeftAssemblyLength=$F[11];RightAssemblyLength=$F[12];AveragePercentIdentityAtLeftPos=$F[13];PercentOfNonUniqueMappingReadsAtLeftPos=$F[14];AveragePercentIdentityAtRightPos=$F[15];PercentOfNonUniqueMappingReadsAtRightPos=$F[16];StartPositionOfConsensusMappingToGenome=$F[17];StartingChromosomeOfConsensusMapping=$F[18];PositionOfTheGenomicMappingOfConsensusStartingPosition=$F[19];EndPositionOfConsensusMappingToGenome=$F[20];EndingChromsomeOfConsnesusMapping=$F[21];PositionOfGenomicMappingOfConsensusEndingPosiiton=$F[22];ConsensusSequences=$F[23];TCHR=na;TSTART=na;SVID=$sv_id;SVType=breakpoint\n";
	}
}
###################################################################




sub usage
{
	die `pod2text $0`;
}
