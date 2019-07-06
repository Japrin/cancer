#!/usr/bin/env perl
#============================================================================
# Name        		: add_mutPos_kmer.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sun Mar 22 22:23:59 2015
# Last Modified By	: 
# Last Modified On	: Sun Mar 22 22:23:59 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl add_mutPos_kmer.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;


if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
	#LTRLCRWDK       0.443   415     WB      LTRLCQWDK       0.310   1747    .       HLA-A31:01      gbm	...	R296Q
	my ($n_peptide,$t_peptide,$aaChange)=@F[0,4,-1];
	my ($a1,$a2)=$aaChange=~/(\D+)\d+(\D+)/;
	my $i=0;
	my $mutPos;
	my $iPos=-1;
	my $jPos=-1;
	my @iAA=();
	my @jAA=();
	my $k=length($n_peptide);
	my @n_peptide=split //,$n_peptide;
	my @t_peptide=split //,$t_peptide;
	for (my $i=0;$i<$k;$i++)
	{
		my $n_aa=$n_peptide[$i];
		my $t_aa=$t_peptide[$i];
		if($n_aa ne $t_aa)
		{
			if($iPos==-1){$iPos=$i;$jPos=$i;}
			else { $jPos=$i; }
			
		}else
		{
			$n_peptide[$i]="\l$n_peptide[$i]";
			$t_peptide[$i]="\l$t_peptide[$i]";
		}
	}
	if(join("",@n_peptide[$iPos..$jPos]) ne $a1 || join("",@t_peptide[$iPos..$jPos]) ne $a2)
	{
		print STDERR "Inconsistant\tmutPos:$iPos..$jPos\t$line\n";
	}
	if($iPos==$jPos) { $mutPos=$iPos+1; }
	else { $mutPos=($iPos+1)."..".($jPos+1); }
	$F[0]=join("",@n_peptide);
	$F[4]=join("",@t_peptide);
	splice (@F,8,0,$mutPos,$k);
	print join("\t",@F)."\n";

}
###################################################################




sub usage
{
    die `pod2text $0`;
}
