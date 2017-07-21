#!/usr/bin/env perl
#============================================================================
# Name        		: TCRasm.compare.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Sat Jun  3 09:47:03 2017
# Last Modified By	: 
# Last Modified On	: Sat Jun  3 09:47:03 2017
# Copyright   		: Copyright (C) 2017
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl TCRasm.compare.pl [option] <tool:infile> <tool:infile1,infile2> ...

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Dumper;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
#my $infile=shift @ARGV;

my %res=();
# tools
my @tnames=();
# sample names
my %snames=();
my $j=0;
foreach (@ARGV){
    $j++;
    my ($tool_name,$ifile)=split /:/,$_;
    my $t_res=processRes($tool_name,$ifile);
    $res{$tool_name}=$t_res;
    push @tnames,$tool_name;
    foreach (keys %{$t_res}){ 
        ## compare the samples in the first input
        if($j==1){$snames{$_}++;}
    }
}
#print Dumper(\%res);
#print Dumper(\%snames);
my $nsample=scalar keys %snames;
my @nDup=(0)x@tnames;
my @nMAIT=(0)x@tnames;
my @nCall=(0)x@tnames;
my %dComp=();
foreach my $sid (sort keys %snames){
    print "#DETAIL\t$sid";
    for(my $iTool=0;$iTool<@tnames;$iTool++){
        my $pVDJ=$res{$tnames[$iTool]}->{$sid}->{"VDJ"};
        my $pCDR3_nt=$res{$tnames[$iTool]}->{$sid}->{"CDR3_nt"};
        my $pCDR3_aa=$res{$tnames[$iTool]}->{$sid}->{"CDR3_aa"};
        for(my $i=0;$i<4;$i++){
            #### ignore difference in allele level
            $pVDJ->[$i]=~s/\*\d+//g;
            ### remove _ and /
            $pVDJ->[$i]=~s/[_\/]//g;
            ### ignore TRBD
            $pVDJ->[$i]=~s/\|.+\|/\|\.\|/;
        }
        ### check duplicates (due to sequencing error; mixcr suffer this)
        if($pVDJ->[0] eq $pVDJ->[1] && $pVDJ->[0] ne "NA"){ 
            $nDup[$iTool]++; $pVDJ->[1]="NA"; 
            #if($iTool==1){ print STDERR "###TEST_A\t$sid\n"; }
        }
        if($pVDJ->[2] eq $pVDJ->[3] && $pVDJ->[2] ne "NA"){ 
            $nDup[$iTool]++; $pVDJ->[3]="NA"; 
            #if($iTool==1){ print STDERR "###TEST_B\t$sid\n"; }
        }
        ### check MAIT
        if($pVDJ->[0]=~/^TRAV1-2\|.+TRAJ(33|20|12)$/ || $pVDJ->[1]=~/^TRAV1-2\|.+TRAJ(33|20|12)$/){ 
            $nMAIT[$iTool]++; 
            #if($iTool==0){ print STDERR "$sid\n"; }
        }
        ### check consistence
        for(my $i=0;$i<4;$i++){
            if($pVDJ->[$i] ne "NA"){
                my $cid=sprintf("%s.%s",$sid,$pVDJ->[$i]);
                $dComp{$cid}+=(2**$iTool);
                $nCall[$iTool]++;
            }
        }
        ### output
        print "\t".join("\t",@$pVDJ);
        #print "\t".join("\t",@$pCDR3_nt);
        #print "\t".join("\t",@$pCDR3_aa);
    }
    print "\n";
    ###print Dumper($res{$tnames[0]})
}
my @v=values(%dComp);
my %compStat=();
foreach (@v){ $compStat{$_}++;  }
print "#COMP\ttoolName\t".join("\t",@tnames)."\n";
print "#COMP\tDupRate\t".join("\t",@nDup)."\t$nsample\n";
print "#COMP\tnMAIT\t".join("\t",@nMAIT)."\t$nsample\n";
print "#VENN\tnCall\t".join("\t",@nCall)."\t".(scalar @v)."\n";
foreach (sort { $a<=>$b } keys %compStat){
    print "#VENN\t$_\t$compStat{$_}\n";
}
foreach (sort keys %dComp){
    print "#DUMP\t$_\t$dComp{$_}\n";
}
#print Dumper(\%dComp);
    
sub processRes
{
    my ($tool_name,$infile)=@_;
    if($tool_name eq "tracer"){
        open $in,$infile or die "Cann't open file $infile ($!) \n";
        my %ret=();
        my $hh=<$in>;
        while(<$in>)
        {
            chomp;
            my $line=$_;
            if(/^\s*$/ || /^#/) { next; }
            my @F=split /\t/;
            my $sid=$F[0];
            my @tID=@F[3,11,19,27];
            my @CDR3_aa=@F[4,12,20,28];
            my @CDR3_nt=@F[5,13,21,29];
            my @VDJ=@F[6,14,22,30];
            my @TPM=@F[7,15,23,31];
            my @productive=@F[8,16,24,32];
            @CDR3_aa=map{ uc } @CDR3_aa;
            @CDR3_nt=map{ uc } @CDR3_nt;
            for (my $i=0; $i <4; $i++){
                if($CDR3_nt[$i] eq "" || $CDR3_nt[$i]=~/Couldn/){ $CDR3_nt[$i]="NA"; }
                if($CDR3_aa[$i] eq "" || $CDR3_aa[$i]=~/Couldn/){ $CDR3_aa[$i]="NA"; }
                #$VDJ[$i]=(split /,/,$VDJ[$i])[0];
                my @tmp=split /\|/,$VDJ[$i];
                $VDJ[$i]=join("|",map { (split /,/,$_)[0] } @tmp);
                $VDJ[$i]=~s/\|N\/A\|/\|\.\|/;
            }
            $ret{$sid}={"CDR3_aa"=>\@CDR3_aa,"CDR3_nt"=>\@CDR3_nt,"VDJ"=>\@VDJ};
        }
        return \%ret;
    }elsif($tool_name eq "mixcr"){
        open $in,$infile or die "Cann't open file $infile ($!) \n";
        my %ret=();
        while(<$in>)
        {
            chomp;
            my $line=$_;
            if(/^\s*$/ || /^#/) { next; }
            my @F=split /\t/;
            my ($pid,$sid,$sfile)=@F[0..2];
            my $sin;
            open $sin,$sfile or die "Cann't open file $sfile ($!) \n";
            my $shh=<$sin>;
            my $n_alpha=0;
            my $n_beta=0;
            my @VDJ_A=();
            my @VDJ_B=();
            my @CDR3_nt_A=();
            my @CDR3_nt_B=();
            my @CDR3_aa_A=();
            my @CDR3_aa_B=();
            while(<$sin>){
                chomp;
                my $sline=$_;
                if(/^\s*$/ || /^#/) { next; }
                my @sF=split /\t/;
                if($sF[8]=~/^TRAC/){
                    $n_alpha++;
                    if($n_alpha<=2){ 
                        ## [$vdj,$cdr3_nt,$cdr3_aa]
                        my $pL=processMiTCRLine(\@sF);
                        push @VDJ_A,$pL->[0];
                        push @CDR3_nt_A,$pL->[1];
                        push @CDR3_aa_A,$pL->[2];
                    }
                }elsif($sF[8]=~/^TRBC/){
                    $n_beta++;
                    if($n_beta<=2){ 
                        my $pL=processMiTCRLine(\@sF);
                        push @VDJ_B,$pL->[0];
                        push @CDR3_nt_B,$pL->[1];
                        push @CDR3_aa_B,$pL->[2];
                    }
                }
            }
            # VDJ
            my @VDJ=();
            while(@VDJ_A<2){ push @VDJ_A,"NA"; }
            while(@VDJ_B<2){ push @VDJ_B,"NA"; }
            push @VDJ,@VDJ_A;
            push @VDJ,@VDJ_B;
            # CDR3_nt
            my @CDR3_nt=();
            while(@CDR3_nt_A<2){ push @CDR3_nt_A,"NA"; }
            while(@CDR3_nt_B<2){ push @CDR3_nt_B,"NA"; }
            push @CDR3_nt,@CDR3_nt_A;
            push @CDR3_nt,@CDR3_nt_B;
            # CDR3_aa
            my @CDR3_aa=();
            while(@CDR3_aa_A<2){ push @CDR3_aa_A,"NA"; }
            while(@CDR3_aa_B<2){ push @CDR3_aa_B,"NA"; }
            push @CDR3_aa,@CDR3_aa_A;
            push @CDR3_aa,@CDR3_aa_B;
            # return data
            $ret{$sid}={"CDR3_aa"=>\@CDR3_aa,"CDR3_nt"=>\@CDR3_nt,"VDJ"=>\@VDJ};
        }
        return \%ret;
    }elsif($tool_name eq "vdjpuzzle"){
        open $in,$infile or die "Cann't open file $infile ($!) \n";
        my %ret=();
        while(<$in>)
        {
            chomp;
            my $line=$_;
            if(/^\s*$/ || /^#/) { next; }
            my @F=split /\t/;
            my ($pid,$sid,$tra_file,$trb_file)=@F;
            my @VDJ=();
            my @CDR3_aa=();
            my @CDR3_nt=();
            my $pTRA=processVDJPuzzleFile($tra_file);
            my $pTRB=processVDJPuzzleFile($trb_file);
            push @CDR3_aa,@{$pTRA->[0]};
            push @CDR3_nt,@{$pTRA->[1]};
            push @VDJ,@{$pTRA->[2]};
            push @CDR3_aa,@{$pTRB->[0]};
            push @CDR3_nt,@{$pTRB->[1]};
            push @VDJ,@{$pTRB->[2]};
            $ret{$sid}={"CDR3_aa"=>\@CDR3_aa,"CDR3_nt"=>\@CDR3_nt,"VDJ"=>\@VDJ};
        }
        return \%ret;
    }
}
###################################################################

sub processVDJPuzzleFile
{
    my ($infile)=@_;
    my $in;
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    my $hh=<$in>;
    my @VDJ=();
    my @CDR3_aa=();
    my @CDR3_nt=();
    my $n=0;
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
        my ($_v,$_d,$_j,$_cdr3nt,$_cdr3aa)=@F[3..7];
        my $_VDJ=join("|",$_v,$_d,$_j);
        $n++;
        if($n<=2){
            push @VDJ,$_VDJ;
            push @CDR3_aa,$_cdr3aa;
            push @CDR3_nt,$_cdr3nt;
        }
    }
    while(@VDJ<2){ push @VDJ,"NA"; }
    while(@CDR3_aa<2){ push @CDR3_aa,"NA"; }
    while(@CDR3_nt<2){ push @CDR3_nt,"NA"; }
    return([\@CDR3_aa,\@CDR3_nt,\@VDJ]);
}

sub processMiTCRLine
{
    my ($pF)=@_;
    my $v=(split /,/,$pF->[5])[0];
    $v=~s/\(.+\)//;
    my $d=(split /,/,$pF->[6])[0];
    if(defined($d)){
        $d=~s/\(.+\)//;
    }
    if(!defined($d) || $d eq "") { $d="."; }
    my $j=(split /,/,$pF->[7])[0];
    $j=~s/\(.+\)//;
    my $vdj=join("|",$v,$d,$j);
    my ($cdr3_nt,$cdr3_aa)=@{$pF}[23,32];
    my $ret=[$vdj,$cdr3_nt,$cdr3_aa];
    return($ret);
}

sub usage
{
    die `pod2text $0`;
}
