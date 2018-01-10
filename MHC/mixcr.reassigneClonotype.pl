#!/usr/bin/env perl
#============================================================================
# Name        		: mixcr.reassigneClonotype.pl
# Author      		: zhenglt@gmail.com
# Version     		: v1.00
# Created On  		: Thu Dec 10 09:10:04 2015
# Last Modified By	: 
# Last Modified On	: Thu Dec 10 09:10:04 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl mixcr.reassigneClonotype.pl [option] <infile> <outfile>

    -s  sampleID [default: SAMPLE]
    -r  reorder by "main alpha > secondary alpha > main beta > secondary beta" [default: OFF]
    -h  display this help and exit

=cut

###-a  TPM threshold of alpha [default: 10]
###-b  TPM threshold of beta [default: 15]

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Table;

my ($in,$out);
my ($opt_h,$opt_s,$opt_r,$opt_a,$opt_b);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s,"r"=>\$opt_r,"a=i"=>\$opt_a,"b=i"=>\$opt_b);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
##if(!defined($opt_a)) { $opt_a=10; }
##if(!defined($opt_b)) { $opt_b=15; }

my $infile=shift @ARGV;
my $outfile=shift @ARGV;

run($infile,$outfile);

###################################################################



sub run
{
    my ($in,$out);
    my ($infile,$outfile)=@_;
    my %clonotype=();
    my %clonotypeCount=();
    my @inputData=();
    my %cellClonotype=();
	my $in_table = Data::Table::fromTSV($infile, 1);
	my $next = $in_table->iterator();
    open $out,">",$outfile or die "Cann't open file $outfile ($!) \n";
	print $out join("\t",$in_table->header,"C_strict","C_strict_b")."\n";
    #chomp $header;
    #if(!$opt_r) { 
    #    ##print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\tmainAlpha\tmainBeta\n"; 
    #    print $out "$header\tFoundMajor\tC_strict\tC_share\tC_strict_b\tC_share_b\tmainAlpha\tmainBeta\n";
    #}
    #else 
    #{ 
    #    print $out "$header\tC_strict\tC_share\tC_strict_b\tC_share_b\n"; 
    #    #print $out "$header\tFoundMajor\tC_strict\tC_share\tC_strict_b\tC_share_b\n";
    #}
    my $cid=0;
    ###while(<$in>)
	while (my $row = $next->())
    {
        #chomp;
        #my $line=$_;
        #if(/^\s*$/ || /^#/) { next; }
        #my @F=split /\t/;
        my $cellID=$row->{"sample"};
        my @VDJ_A=($row->{"VDJ_A1"},$row->{"VDJ_A2"});
        my @VDJ_B=($row->{"VDJ_B1"},$row->{"VDJ_B2"});
        my @CDR3_aa_A=($row->{"CDR3_aa_A1"},$row->{"CDR3_aa_A2"});
        my @CDR3_aa_B=($row->{"CDR3_aa_B1"},$row->{"CDR3_aa_B2"});
        my @count_A=($row->{"count_A1"},$row->{"count_A2"});
        my @count_B=($row->{"count_B1"},$row->{"count_B2"});

        my $mainA=-1;
        my $mainB=-1;
        
        #### main alpha
        if($CDR3_aa_A[0] ne "NA" && ($VDJ_A[1] eq "NA" || $CDR3_aa_A[1] eq "NA")){
            $mainA=0;
        }elsif($CDR3_aa_A[1] ne "NA" && ($VDJ_A[0] eq "NA" || $CDR3_aa_A[0] eq "NA")){
            $mainA=1;
        }elsif($CDR3_aa_A[0] ne "NA" && $CDR3_aa_A[1] ne "NA"){
            $mainA=$count_A[0]>$count_A[1]?0:1;
        }
        #### main beta
        if($CDR3_aa_B[0] ne "NA" && ($VDJ_B[1] eq "NA" || $CDR3_aa_B[1] eq "NA")){
            $mainB=0;
        }elsif($CDR3_aa_B[1] ne "NA" && ($VDJ_B[0] eq "NA" || $CDR3_aa_B[0] eq "NA")){
            $mainB=1;
        }elsif($CDR3_aa_B[0] ne "NA" && $CDR3_aa_B[1] ne "NA"){
            $mainB=$count_B[0]>$count_B[1]?0:1;
        }
        my $foundMajor=1;
        if($mainA==-1 || $mainB==-1)
        {
            print STDERR "No major: mainAlpha=$mainA\tmainBeta=$mainB\t$cellID\n";
            $mainA=0;
            $mainB=2;
            $foundMajor=0;
            #next;
        }
        #### if one productive alpha and one productive beta is the same, belong to the same clone
        my @ID_A_forClone=();
        my @ID_B_forClone=();
        if($CDR3_aa_A[0] ne "NA") { push @VDJ_A,$CDR3_aa_A[0]; }
        if($CDR3_aa_A[1] ne "NA") { push @VDJ_A,$CDR3_aa_A[1]; }
        if($CDR3_aa_B[0] ne "NA") { push @VDJ_B,$CDR3_aa_B[0]; }
        if($CDR3_aa_B[1] ne "NA") { push @VDJ_B,$CDR3_aa_B[1]; }
        my $beCounted=0;
        my $beFound=0;
        foreach my $_id_A (@ID_A_forClone){
            foreach my $_id_B (@ID_B_forClone){
                if($clonotype{"$_id_A:$_id_B"}){
                    ###$cellClonotype{$cellID}={'nid'=>$clonotype{"$_id_A:$_id_B"}};
                    if($beCounted==0){
                        ### move here
                        $cellClonotype{$cellID}={'nid'=>$clonotype{"$_id_A:$_id_B"}};
                        $clonotypeCount{$cellClonotype{$cellID}->{'nid'}}++;
                        $beCounted=1;
                    }
                    $beFound=1;
                }else{
                }
            }
        }
        if(!$beFound){
            $cid++;
            my $nid=sprintf("C%06d",$cid);
            if(@ID_A_forClone>0 && @ID_B_forClone>0){
                $cellClonotype{$cellID}={'nid'=>$nid};
                $clonotypeCount{$nid}=1;
                foreach my $_id_A (@ID_A_forClone){
                    foreach my $_id_B (@ID_B_forClone){
                        $clonotype{"$_id_A:$_id_B"}=$nid;
                    }
                }
            }else{
                $cellClonotype{$cellID}={'nid'=>"NoPPair"};
                #$clonotypeCount{$nid}=1;
            }
        }
        #### 
    }
	$next = $in_table->iterator();
	for(my $i = 0; $i < $in_table->nofRow; $i++)
    {
		my @F=$in_table->row($i);
        my $cellID=$F[0];
        my $nid="NA";
        my $strMainAlpha="NA";
        my $strMainBeta="NA";
        my $nid_count=0;
        $nid=$cellClonotype{$cellID}->{'nid'};
        #$strMainAlpha=$cellClonotype{$cellID}->{'mainAlpha'};
        #$strMainBeta=$cellClonotype{$cellID}->{'mainBeta'};
        $nid_count=$clonotypeCount{$nid};
        if(!defined($nid_count)){ $nid_count=0; }
        { 
            print $out join("\t", @F,"${opt_s}_$nid:$nid_count",$nid_count>1?"Clonal":"NoClonal")."\n"; 
        }
    }
}

sub usage
{
    die `pod2text $0`;
}
