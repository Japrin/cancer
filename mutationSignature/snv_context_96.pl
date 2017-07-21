#!/usr/bin/env perl
#============================================================================
# Name        		: snv_context_96.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Tue Jun 24 22:53:31 2014
# Last Modified By	: 
# Last Modified On	: Tue Jun 24 22:53:31 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl snv_context_96.pl [option] [infile]

	-a	chr column [0-based, default 0] 
	-b	pos column [0-based, default 1] 
	-i	ref column [0-based, default 3] 
	-j	mut column [0-based, default 4] 
	-f	reference sequence in fasta format [default: /DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12/human_g1k_v37_decoy.fasta ]
	-s	sampleID [default: SAMPLE]
	-m	0-based sample column index [default: not used, assum infile contain only one sample]
	-h	display this help and exit
    Ignore the lines start with "#"

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Bio::DB::Sam;
use File::Basename;

my ($in,$out);
my ($opt_h,$opt_a,$opt_b,$opt_i,$opt_j,$opt_f,$opt_s,$opt_m);
GetOptions("h"	=>\$opt_h,"a=i"=>\$opt_a,"b=i"=>\$opt_b,"i=i"=>\$opt_i,"j=i"=>\$opt_j,"f=s"=>\$opt_f,"s=s"=>\$opt_s,"m=i"=>\$opt_m);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_a)) { $opt_a=0; }
if(!defined($opt_b)) { $opt_b=1; }
if(!defined($opt_i)) { $opt_i=3; }
if(!defined($opt_j)) { $opt_j=4; }
if(!defined($opt_f)) { $opt_f="/DBS/DB_temp/zhangLab/broad/bundle/2.8/b37/bwa_0.7.12/human_g1k_v37_decoy.fasta"; }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
if(!defined($opt_m)) { $opt_m=-1; }
#my $infile=shift @ARGV;
#my $outfile=shift @ARGV;

my %context=();
my $REF=readReference($opt_f);
##my $fai = Bio::DB::Sam::Fai->load($opt_f);
my $iRecord=0;
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#print STDERR "$opt_i\t$opt_j\t$F[$opt_i]\t$F[$opt_j]\t----\n";
	my ($chr,$pos,$from,$to)=@F[$opt_a,$opt_b,$opt_i,$opt_j];
	$to=(split /,/,$to)[0];
    my $sampleID=$opt_s;
    if($opt_m != -1){ $sampleID=$F[$opt_m]; }
	if(length($from) == length($to) && length($from) == 1)
	{
		#snv
        if(!$REF->{$chr}){ print STDERR "ERROR: $chr not in reference ($line)\n"; next; }
        my $chrLen=length($REF->{$chr});
        if($pos<=1 || $pos >=$chrLen){ print STDERR "ERROR: $pos out of range of $chr (1:$chrLen) ($line)\n"; next; }
        my $dna_string=substr($REF->{$chr},$pos-2,3);
        #my $dna_string = $fai->fetch(sprintf("$chr:%d-%d",$pos-1,$pos+1));
		$dna_string="\U$dna_string";
		## validation
		my @ss=split //,$dna_string;
		if($ss[1] ne $from)
		{
			print STDERR "ERROR: $from no equal to reference ($line)\n";
			next;
		}
		if($from ne "T" && $from ne "C")
		{
			$from=~tr/ATCG/TAGC/;
			$to=~tr/ATCG/TAGC/;
			$dna_string=~tr/ATCG/TAGC/;
			$dna_string=reverse($dna_string);
			@ss=split //,$dna_string;
		}
        ###$context{sprintf("%s(%s>%s)%s",$ss[0],$from,$to,$ss[2])}++;
		$context{$sampleID}->{sprintf("%s(%s>%s)%s",$ss[0],$from,$to,$ss[2])}++;
        $iRecord++;
        if($iRecord%1000==0){
            printf STDERR "processed record: %d\n",$iRecord;
        }
	}
}
my @allSamples=sort(keys(%context));
print "context\t".join("\t",@allSamples)."\n";
#foreach (sort keys %context)
foreach my $mut ("C>A","C>G","C>T","T>A","T>C","T>G")
{
	foreach my $f5 ("A","C","G","T")
	{
		foreach my $f3 ("A","C","G","T")
		{
			my $k="$f5($mut)$f3";
            print "$k";
            foreach (@allSamples){
                print "\t".($context{$_}->{$k}?$context{$_}->{$k}:0);
            }
            print "\n";
            ##print "$k\t".($context{$k}?$context{$k}:0)."\n";
		}
	}
}
###################################################################

sub readReference
{
    my ($infile)=@_;
    my %R=();
    my $in;
    if($infile=~/\.gz$/){ open $in, "gzip -cd $infile |" or die "$!"; }
    else{ open $in, $infile or die "$!"; }
    my $seq_name="";
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        if($line=~/^>/){
            ($seq_name)=$line=~/^>(.+?)(\s|$)/;
        }elsif($seq_name ne "") {
            $R{$seq_name}.=$line;
        }
    }
    return \%R;
}


sub usage
{
	die `pod2text $0`;
}
