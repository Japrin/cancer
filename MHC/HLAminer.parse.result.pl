#!/usr/bin/perl
#============================================================================
# Name        		: HLAminer.parse.result.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Tue Dec 30 08:27:18 2014
# Last Modified By	: 
# Last Modified On	: Tue Dec 30 08:27:18 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl HLAminer.parse.result.pl [option] <infile>
    -v	verbose output [default OFF]
    -s	sampleID [default "sampleID" ]

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_v,$opt_s);
GetOptions("h"	=>\$opt_h,"v"=>\$opt_v,"s=s"=>\$opt_s);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_s)) { $opt_s="sampleID"; }


####HLA-A
####	Prediction #1 - A*03
####		A*03:02P,1196980.848,9.43e-06,50.3
####	Prediction #2 - A*24
####		A*24:54,740899.6628,7.37e-06,51.3
####		A*24:02P,740899.6628,7.37e-06,51.3
####		A*24:61,419226.8072,7.37e-06,51.3
####	Prediction #3 - A*11
####		A*11:01P,553244.384,1.89e-05,47.2
####
####HLA-B
####	Prediction #1 - B*47
####		B*47:01P,469979.637,2.76e-12,115.6
####	Prediction #2 - B*27
####		B*27:96P,2250,7.52e-05,41.2
####		B*27:97,2250,7.52e-05,41.2
####		B*27:99,2250,7.52e-05,41.2
####		B*27:82,2250,7.52e-05,41.2
####		B*27:98,2250,7.52e-05,41.2
####		B*27:84,2250,7.52e-05,41.2
####		B*27:90P,2250,7.52e-05,41.2
####
####HLA-DPB1
####	Prediction #1 - DPB1*23
####		DPB1*23:01P,7200,5.89e-04,32.3
####	Prediction #1 - DPB1*04 (same score as above)
####		DPB1*04:01P,7200,5.89e-04,32.3
####		DPB1*04:02P,7200,5.89e-04,32.3
####	Prediction #1 - DPB1*39 (same score as above)
####		DPB1*39:01,7200,5.89e-04,32.3
####	Prediction #1 - DPB1*121 (same score as above)
####		DPB1*121:01,7200,5.89e-04,32.3
####
####HLA-DQA1

my %res=();
open $in,$infile or die "Cann't open file $infile ($!) \n";
if($opt_v)
{
    print STDERR "#### HLAminer output ####\n";
    print STDERR "Output\tGene\tPred_#\tAllele_Group\tProtein_Coding_Group\tScore\tEval\tConfidence\n";
}
while(<$in>)
{
    chomp;
    my $line=$_;
    my $gene="";
    if(/^HLA/)
    {
    	$gene=$_;
	while(<$in>)
	{
		chomp;
		my $isCont=1;
		if(/^\tPrediction #(\d+) - (.+)$/)
		{
			my $pred_i=$1;
			my $alleleGroup=$2;
			$alleleGroup=~s/ \(same score as above\)//;
			while(<$in>)
			{
				chomp;
				if(/^\t\t(.+)$/)
				{
					my $proteinCodingGroupLine=$1;
					my ($proteinCodingGroup,$score,$eval,$confidence)=split /,/,$proteinCodingGroupLine;
					if($opt_v)
					{
						print STDERR join("\t",("HLAminer",$gene,$pred_i,$alleleGroup,$proteinCodingGroup,$score,$eval,$confidence))."\n";
					}
					if($pred_i<3)
					{
						if(exists($res{$gene}->{$pred_i}))
						{
							#set ambiguity
							my @k=keys %{$res{$gene}->{$pred_i}};
							my @kk=keys %{$res{$gene}->{$pred_i}->{$k[0]}};
							if($k[0] ne $alleleGroup) 
							{
								if($res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'score'} == $score)
								{
									$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'ambiguity'}="2digit_ambiguity";
								}
							}else
							{
								if($res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'score'} == $score)
								{
									$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'ambiguity'}="4digit_ambiguity";
								}

							}
						}else
						{
							
							$res{$gene}->{$pred_i}->{$alleleGroup}->{$proteinCodingGroup}={'score'=>$score,'eval'=>$eval,'confidence'=>$confidence,'ambiguity'=>"NoAmbiguity"};
						}
					}


				}elsif(/^\tPrediction #(\d+) - (.+)$/)
				{
					$pred_i=$1;
					$alleleGroup=$2;
					$alleleGroup=~s/ \(same score as above\)//;
				}else
				{
					$isCont=0;
					last
				}
			}

		}else
		{
			$isCont=0;
		}
		if(!$isCont)
		{
			last;
		}

	}

    }
}
print STDERR "#### filtering output ####\n";
print STDERR "Output\tGene\tPred_#\tAllele_Group\tProtein_Coding_Group\tScore\tEval\tConfidence\tAmbiguity\n";
#foreach my $gene (sort keys %res)
print "sampleID\tHLA-A\tHLA-B\tHLA-C\tHLA-E\tHLA-F\tHLA-G\tHLA-DPA1\tHLA-DPB1\tHLA-DQA1\tHLA-DQB1\tHLA-DRA\tHLA-DRB1\tHLA-DRB2\tHLA-DRB3\tHLA-DRB4\tHLA-DRB5\tHLA-DRB6\tHLA-DRB7\tHLA-DRB8\tHLA-DRB9\n";
print "$opt_s";
foreach my $gene ("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB2","HLA-DRB3","HLA-DRB4","HLA-DRB5","HLA-DRB6","HLA-DRB7","HLA-DRB8","HLA-DRB9")
{
	my @gt=();
	foreach my $pred_i (keys %{$res{$gene}})
	{
		my @k=keys %{$res{$gene}->{$pred_i}};
		my @kk=keys %{$res{$gene}->{$pred_i}->{$k[0]}};
		if($opt_v)
		{
			print STDERR "Filtering\t$gene\t$pred_i\t$k[0]\t$kk[0]\t".join("\t",$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'score'},
								$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'eval'},
								$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'confidence'},
								$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'ambiguity'})."\n";
		}
		my $ab=$res{$gene}->{$pred_i}->{$k[0]}->{$kk[0]}->{'ambiguity'};
		if($ab eq "NoAmbiguity")
		{
			push @gt,$kk[0];
		}elsif($ab eq "4digit_ambiguity")
		{
			push @gt,$k[0];
		}elsif($ab eq "2digit_ambiguity")
		{
			push @gt,"NR";
		}
		
		

	}
	print "\t".join(",",@gt);

}
print "\n";
###################################################################




sub usage
{
    die `pod2text $0`;
}


