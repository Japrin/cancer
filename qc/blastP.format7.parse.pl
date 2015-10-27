#!/usr/bin/env perl
#============================================================================
# Name        		: blastP.format7.parse.pl
# Author      		: zeminz_pkuhpc
# Version     		: v1.00
# Created On  		: Fri Jul  3 19:05:34 2015
# Last Modified By	: 
# Last Modified On	: Fri Jul  3 19:05:34 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl blastP.format7.parse.pl [option] <infile>

    -v  verbose output
    -s  sample id [dedault SAMPLE]
    -c  min coverage [dedault 0.90]
    -p  min identity [dedault 0.95]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_v,$opt_s,$opt_c,$opt_p);
GetOptions("h"	=>\$opt_h,"v"=>\$opt_v,"s=s"=>\$opt_s,"c=f"=>\$opt_c,"p=f"=>\$opt_p);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
if(!defined($opt_c)) { $opt_c=0.9; }
if(!defined($opt_p)) { $opt_p=0.95; }


if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
$/="hits found";
my $hhh=<$in>;
my $n_total=0;
my $n_hit=0;
while(<$in>)
{
    #chomp;
    #$/="\n";
    $n_total++;
    my @lines=split /\n/;
    my $n=@lines;
    for(my $i=0; $i<$n; $i++)
    {
        my @F=split /\t/,$lines[$i];
        if(@F<17) { next; }
        ### the top1
        my ($query_id,$subject_id,$subject_tax_id,$subject_sci_name,$subject_com_name,$query_len,$subject_len,$p_identity,$alignment_len)=@F[0..8];
        if($alignment_len/$query_len>$opt_c && $p_identity>$opt_p)
        {
            $n_hit++;
            if($opt_v) 
            {
                print "#$opt_s\t".join("\t",@F)."\n";
            }
        }
        #print join("\t",@F)."\n";
        last;
    }
    #1  query id
    #2  subject id
    #3  subject tax ids
    #4   subject sci names
    #5  subject com names
    #6  query length
    #7  subject length
    #8  %identity
    #9  alignment length
    #10 mismatches
    #11 gap opens
    #12 q. start
    #13 q. end
    #14 s. start
    #15 s. end
    #16 evalue
    #17 bit score
    # Fields: query id, subject id, subject tax ids, subject sci names, subject com names, query length, subject length, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    # # 23 hits found
    # HWI-ST1352:307:HV2V2ADXX:2:1101:8595:2233       gi|685048051|emb|LN596563.1|    7962    Cyprinus carpio common carp     53      179080  100.00  36      0       0       1       36      141655  141690  8e-09   67.6

}
print "#Sample\ttotal\thit\thit rate\n";
printf "$opt_s\t$n_total\t$n_hit\t%4.4f\n",$n_hit/$n_total;
###################################################################




sub usage
{
    die `pod2text $0`;
}
