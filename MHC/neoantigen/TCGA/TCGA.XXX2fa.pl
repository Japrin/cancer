#!/usr/bin/perl
#============================================================================
# Name        		: XXX2fa.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Mon Oct 13 10:41:57 2014
# Last Modified By	: 
# Last Modified On	: Mon Oct 13 10:41:57 2014
# Copyright   		: Copyright (C) 2014
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl XXX2fa.pl [option] <infile>
    -e  ensemble fa [default: "/DBS/DB_temp/zhangLab/ensemble/release69/homo_sapiens/fasta/pep/Homo_sapiens.GRCh37.69.pep.all.fa.gz" ]
    -o  output prefix [default: ./out ]
    -l  pep length [default: 9]
    -i  column number (0-based) of chr,pos,hugo_symbol,enst_id,aaChang,varClass,sampleID,ref,alt in XXX file [default: 4,5,0,39,47,8,15,35,36 ]	
    -m  max output sequence length [default 20000]
    -s  output sequence name's prefix  [default "sample"]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path qw(make_path remove_tree);

my ($in,$out,$out1,$out2,$out3,$out4);
my ($opt_h,$opt_e,$opt_o,$opt_l,$opt_m,$opt_s,$opt_i);
GetOptions("h"	=>\$opt_h,"e=s"=>\$opt_e,"o=s"=>\$opt_o,"l=i"=>\$opt_l,"m=i"=>\$opt_m,"s=s"=>\$opt_s,"i=s"=>\$opt_i);
if(@ARGV<1 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_e)) { $opt_e="/DBS/DB_temp/zhangLab/ensemble/release69/homo_sapiens/fasta/pep/Homo_sapiens.GRCh37.69.pep.all.fa.gz"; }
if(!defined($opt_o)) { $opt_o="./out"; }
if(!defined($opt_l)) { $opt_l=9; }
if(!defined($opt_m)) { $opt_m=20000; }
if(!defined($opt_s)) { $opt_s="sample"; }
if(!defined($opt_i)) { $opt_i="4,5,0,39,47,8,15,35,36"; }
my @opt_i=split /,/,$opt_i;

my ($outDir)=$opt_o=~/(.+)\//;
make_path($outDir);

my %proSeq=();
readEnsembleProteinFa(\%proSeq,$opt_e);

open $in,$infile or die "Cann't open file $infile ($!) \n";
open $out1,">","$opt_o.ori.fa" or die "$!";
open $out2,">","$opt_o.mut.fa" or die "$!";
open $out3,">","$opt_o.bed" or die "$!";
open $out4,">","$opt_o.used.XXX" or die "$!";

my $gID=1;
my $gBeg=0;
my $gEnd=0;
my $gSeq1="";
my $gSeq2="";
#my $paddingSeqLen=$opt_l;
#my $paddingSeq="A"x$paddingSeqLen;

while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
	if(/^Hugo_Symbol/) { print $out4 "$line\n"; next; }
    my @F=split /\t/;
    my ($chr,$pos,$hugo_symbol,$enst_id,$aaChange,$varClass,$sampleID,$ref,$alt)=@F[@opt_i];
    if($varClass ne "Missense_Mutation") { next; }
    #print join("\t",($chr,$pos,$hugo_symbol,$enst_id,$aaChange,$varClass))."\n";
    my ($a1,$aPos,$a2)=$aaChange=~/p.(\D+)(\d+)(\D+)$/;     #SNP or DNP or TNP ...
    if(!defined($aPos))
    {
        print STDERR "ERROR\tparse error\t$line\n";
        next;
    }
    my $subLen=length($a1);
    my $o1Seq="";
    my $o2Seq="";
    my $tA="";
	my $mutPosInPeptide=-1;		## 1-based coordinate

    if($proSeq{$enst_id}->{'length'}<$aPos+$subLen-1)
    {
        print STDERR "ERROR\tmaf annotation error\t$line\t$proSeq{$enst_id}->{'length'}\t$proSeq{$enst_id}->{'seq'}\n";
        next;
    }elsif($proSeq{$enst_id}->{'length'}<$opt_l)
    {
        print STDERR "ERROR\tun-availabe prediction\t$line\t$proSeq{$enst_id}->{'length'}\t$proSeq{$enst_id}->{'seq'}\n";
        next;
    }else
    {
        if($proSeq{$enst_id}->{'length'}<$aPos+$subLen-1+($opt_l-1) && $aPos<$opt_l )
        {
            ## short sequence
            $o1Seq=$proSeq{$enst_id}->{'seq'};
            $o2Seq=$o1Seq;
            substr($o2Seq,$aPos-1,$subLen,$a2);
			$mutPosInPeptide=$aPos;
        }elsif($proSeq{$enst_id}->{'length'}<$aPos+$subLen-1+($opt_l-1))
        {
            ## near end of sequence
            $o1Seq=substr($proSeq{$enst_id}->{'seq'},$aPos-($opt_l-1)-1);
            $o2Seq=$o1Seq;
            substr($o2Seq,($opt_l-1),$subLen,$a2);
			$mutPosInPeptide=$opt_l;
        }elsif($aPos<$opt_l)
        {
            ## near beg of sequence
            $o1Seq=substr($proSeq{$enst_id}->{'seq'},0,$aPos+$subLen-1+($opt_l-1));
            $o2Seq=$o1Seq;
            substr($o2Seq,$aPos-1,$subLen,$a2);
			$mutPosInPeptide=$aPos;
        }else
        {
            ## mid part
            $o1Seq=substr($proSeq{$enst_id}->{'seq'},$aPos-($opt_l-1)-1,$opt_l+$subLen-1+($opt_l-1));
            $o2Seq=$o1Seq;
            substr($o2Seq,($opt_l-1),$subLen,$a2);
			$mutPosInPeptide=$opt_l;
        }
        $tA=substr($proSeq{$enst_id}->{'seq'},$aPos-1,$subLen);
        if($tA ne $a1)
        {
            print STDERR "Discordant:\t$line\t$proSeq{$enst_id}->{'length'}\t$proSeq{$enst_id}->{'seq'}\n";
        }
    }
	my ($o1_padding_n1,$o1_padding_n2,$o2_padding_n1,$o2_padding_n2);
	($o1Seq,$o1_padding_n1,$o1_padding_n2)=addPadding($o1Seq,$mutPosInPeptide,$opt_l);
	($o2Seq,$o2_padding_n1,$o2_padding_n2)=addPadding($o2Seq,$mutPosInPeptide,$opt_l);

	my $pipetideLen1=length($o1Seq);
	#my $pipetideLen2=length($o2Seq);
	my $o_name="$sampleID.$chr:$pos.$ref.$alt.$hugo_symbol.$enst_id.$proSeq{$enst_id}->{'ensp_id'}.$aaChange";
	my $o_beg=-1;
	my $o_end=-1;
	if($gEnd+$pipetideLen1>$opt_m)
	{
		print $out1 ">${opt_s}_${gID}\n$gSeq1\n";
		print $out2 ">${opt_s}_${gID}\n$gSeq2\n";
		$gID++;
		$gSeq1=$o1Seq;
		$gSeq2=$o2Seq;
		$o_beg=0;
		$o_end=$o_beg+$pipetideLen1;
		$gEnd=0+$pipetideLen1;

	}else
	{
		$o_beg=$gEnd;
		$o_end=$o_beg+$pipetideLen1;
		$gEnd=$gEnd+$pipetideLen1;
		$gSeq1.=$o1Seq;
		$gSeq2.=$o2Seq;
	}
	#$o_beg=$o_beg+$mutPosInPeptide-1;
	#$o_end=$o_beg+$subLen;
	print $out3 "${opt_s}_$gID\t".($o_beg+$o1_padding_n1)."\t".($o_beg+$opt_l+$opt_l-1-$o1_padding_n2)."\t$o_name\t$o1_padding_n1\t$o1_padding_n2\n";
	
	####print $out1 ">$sampleID.$chr:$pos.$ref.$alt.$hugo_symbol.$enst_id.$proSeq{$enst_id}->{'ensp_id'}.$aaChange\n$o1Seq"."\n";
	####print $out2 ">$sampleID.$chr:$pos.$ref.$alt.$hugo_symbol.$enst_id.$proSeq{$enst_id}->{'ensp_id'}.$aaChange\n$o2Seq"."\n";
	print $out4 "$line\n";
}
if($gSeq1 ne "")
{
	print $out1 ">${opt_s}_${gID}\n$gSeq1\n";
	print $out2 ">${opt_s}_${gID}\n$gSeq2\n";
}
###################################################################

sub addPadding
{
	my ($seq,$mutPos,$peptideLen)=@_;
	###  peptideLen=9
	###  length(seq)=14
	###  mutPos=8
	###  -------X------
	### . n1=1
	###            x=6
	###                ..  n2=2  
	my $n1=$peptideLen-$mutPos;
	my $x=length($seq)-$mutPos;
	my $n2=$peptideLen-1-$x;
	return ("A"x$n1.$seq."A"x$n2,$n1,$n2);
}


sub readEnsembleProteinFa
{
    my $in;
    my ($pList,$infile)=@_;
    if($infile=~/\.gz$/) { open $in,"gzip -cd $infile|" or die "cann't open file $infile ($!) \n"; }
    else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }

    my $ensp_id="";
    my $enst_id="";
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        if(/^>/)
        {
            # a new protein
            #>ENSP00000381008 pep:known chromosome:GRCh37:19:8959608:9091814:-1 gene:ENSG00000181143 transcript:ENST00000397910
            ($ensp_id,$enst_id) = /^>(.+?)\s.+transcript:(.+?)\s/;
            if(defined($ensp_id) && defined($enst_id))
            {
                $pList->{$enst_id}={'length'=>0,'seq'=>"",'ensp_id'=>$ensp_id};
                #print STDERR "OK:\t$ensp_id\t$enst_id\n";
            }else
            {
                $ensp_id="";
                $enst_id="";
                print STDERR "ERROR: $line\n";
            }
        }else
        {
            if(exists($pList->{$enst_id}))
            {
                $pList->{$enst_id}->{'length'} += length($line);
                $pList->{$enst_id}->{'seq'} .= $line;
            }else
            {
                print STDERR "ERROR: $enst_id\t$line\n";
            }
        }

    }
}

sub usage
{
    die `pod2text $0`;
}
