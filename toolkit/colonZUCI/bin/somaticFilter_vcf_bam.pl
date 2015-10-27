#!/usr/bin/perl
#============================================================================
# Name        		: somaticFilter_vcf_bam.pl
# Author      		: japrin
# Version     		: v1.00
# Created On  		: Wed Sep  7 10:14:48 2011
# Last Modified By	: 
# Last Modified On	: Wed Sep  7 10:14:48 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl somaticFilter_vcf_bam.pl [option] <tumor SNP file (vcf)> <normal file (Bam)>

	-f=<file>		index fasta file
	log=<file>		log file
	--max_alt_depth=<int>	max alt bases in normal allowed. > <int> will be filtered [default 2]
	--max_alt_freq=<float>	max alt bases freq in normal allowed. [default 0.1]
	--MQ			min mapping quality [default 10]
	--BQ			min base quality [default 20]
	--snv			filter those in dbSNP or 1KG
	-h			display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
#use lib "/ifshk1/BC_CANCER/01bin/lib/perl/local/lib/perl5/site_perl";
#use lib "/ifshk1/BC_CANCER/01bin/lib/perl/local/lib/perl5/x86_64-linux-thread-multi/"; 
#use Statistics::Descriptive;
#use Text::NSP::Measures::2D::Fisher::twotailed;
#use Statistics::PointEstimation;
#use Statistics::TTest;
use lib "/home/zhoudonger/01bin/lib/perl/lib/perl/5.10.1";
use Bio::DB::Sam;
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_f,$opt_log,$opt_BQ,$opt_MQ,$opt_max_dep,$opt_max_freq,$opt_snv);
GetOptions("h"	=>\$opt_h,"f=s"=>\$opt_f,"log=s"=>\$opt_log,"MQ=i"=>\$opt_MQ,"BQ=i"=>\$opt_BQ,"max_alt_depth=f"=>\$opt_max_dep,"max_alt_freq=f"=>\$opt_max_freq,"snv"=>\$opt_snv);
if(@ARGV<2 || $opt_h) { usage(); }
my $vcffile=shift @ARGV;
my $bamFile=shift @ARGV;

if(!defined($opt_MQ)) { $opt_MQ=10; }
if(!defined($opt_BQ)) { $opt_BQ=20; }
if(!defined($opt_max_dep)) { $opt_max_dep=2; }
if(!defined($opt_max_freq)) { $opt_max_freq=0.1; }
if($opt_log) { open $out,">",$opt_log or die "Cann't open file $opt_log ($!)\n"; }
## read and parse bam file
my ($sam);
if(defined($opt_f) && -e $opt_f) { $sam=Bio::DB::Sam->new(-bam=>$bamFile,-fasta=>$opt_f); }
else { $sam=Bio::DB::Sam->new(-bam=>$bamFile); }

my ($g_chr,$g_pos);
my @g_line=();
my %g_Out=();
my $callback=sub
{
	my ($seqid,$pos,$p) = @_;
	if($seqid ne $g_chr || $pos ne $g_pos) { return; }
	doPileup($seqid,$pos,$p,$sam,\%g_Out,0);
};

if($vcffile=~/\.gz$/) { open $in,"bgzip -cd $vcffile | " or die "Cann't open file $vcffile ($!) \n"; }
else { open $in,$vcffile or die "Cann't open file $vcffile ($!) \n"; }
while(<$in>)
{
	chomp;
	if(/^\s*$/ || /^#/) { print "$_\n"; next; }
	#chr1    26102807        .       G       A       62.15   PASS    AAChange=NM_001145454:c.C98T:p.P33L;AC=1;AF=0.50;AN=2;BaseQRankSum=1.398;Conserved=712(lod=1008);DP=11;Dels=0.00;ExonicFunc=nonsynonymous SNV;FS=5.441;Func=exonic;Gene=STMN1;HRun=0;HaplotypeScore=0.0000;LJB_LRT=0.999912;LJB_MutationTaster=0.966196;LJB_PhyloP=0.999631;MQ=37.00;MQ0=0;MQRankSum=-0.597;PolyPhen2=0.025;QD=5.65;ReadPosRankSum=0.278;SB=-10.63;SIFT=0.01;VQSLOD=4.0341;ensGene=ENSG00000117632:exonic:nonsynonymous SNV:ENST00000236291:c.C98T:p.P33L       GT:AD:DP:GQ:PL  0/1:7,4:11:92.14:92,0,170
	my @field=split /\t/;
	my ($chr,$pos,$rsID,$ref,$alt)=@field[0,1,2,3,4];
	if($rsID ne "." && $opt_snv) { next; }
	if(/1000G.*?=(.+?);/ && $opt_snv) { next; }
	my $beg=$pos;
	my $end=$pos;
	$g_chr=$chr;
	$g_pos=$pos;
	initGOut(\%g_Out,$alt);
	$sam->pileup("$chr:$beg-$end",$callback);
	if($g_Out{'refDepth'}+$g_Out{'ML_ALT_Depth'}==0) { printf STDERR "uncovered in normal:\t$_\n"; next; }
	if(($g_Out{'ML_ALT_Depth'}/($g_Out{'refDepth'}+$g_Out{'ML_ALT_Depth'}))>$opt_max_freq || $g_Out{'ML_ALT_Depth'}>$opt_max_dep)
	{
		printf STDERR "ref:%s\talt_depth:%s\n",$g_Out{'refDepth'},$g_Out{'ML_ALT_Depth'};
	}else
	{
		print "$_\n";
	}
	%g_Out=();
}

############################################################################
sub initGOut
{
	my ($pList,$alt)=@_;
	$pList->{'base4FET'}=[];
	$pList->{'refBaseQual'}=[];
	$pList->{'altBaseQual'}=[];
	$pList->{'refMQ'}=[];
	$pList->{'altMQ'}=[];
	$pList->{'refStrand'}=[];
	$pList->{'altStrand'}=[];
	$pList->{'indel'}=0;
	$pList->{'MQ0'}=0;
	$pList->{'refBase'}="N";
	$pList->{'ML_ALT'}=$alt;
	$pList->{'refDepth'}=0;
	$pList->{'ML_ALT_Depth'}=0;
	$pList->{'bqMean'}=-1;
	$pList->{'mqMean'}=-1;
	$pList->{'refBqMean'}=-1;
	$pList->{'refMqMean'}=-1;
	$pList->{'distMean'}=-1;
	$pList->{'strandBiasP'}=-1;
	$pList->{'baseQualP'}=-1;
	$pList->{'MQP'}=-1;
	$pList->{'A'}=0;
	$pList->{'T'}=0;
	$pList->{'C'}=0;
	$pList->{'G'}=0;
	$pList->{'N'}=0;
}
sub doPileup
{
	my ($seqid,$pos,$p,$sam,$pList,$isTumor) = @_;
	my $refbase = $sam->segment($seqid,$pos,$pos)->dna;
	if($out) { print $out "$seqid\t$pos\t$refbase:\n"; }
	my @_d=();
	my @_base=();
	my @_bq=();
	my @_mq=();
	my @_strand=();
	my $_altDepth=0;
	my $_refDepth=0;
	my $ML_ALT_Dep=0;
	my $ML_ALT="N";
	my ($_dist_mean,$_bq_mean,$_mq_mean)=(0,0,0);
	
	for my $pileup (@$p)
	{
		my $b =$pileup->alignment;
		if($pileup->indel || $pileup->is_refskip)
		{	# don't deal with these ;-)
			$pList->{'indel'}++;
			next;
		}
		my $_read_pos=$pileup->pos;		#not cycle!!
		my $readID=$b->query->name;
		my $qbase=substr($b->qseq,$pileup->qpos,1);
		next if $qbase =~ /[nN]/;
		my $qscore=$b->qscore->[$pileup->qpos];
		my $mq=$b->qual;
		my $strand=$b->strand;
		my $cigar=$b->cigar_str;
		my $q_start=$b->query->start;
		my $q_end=$b->query->end;
		#my $r_length=length($b->qseq);
		#my $q_length=$b->query->length;
		#my @_c=$cigar=~/\*|(?:([0-9]+)([MIDNSHPX=]))/g;
		my $_dist;
		if($strand == 1)
		{
			$_dist=$q_end-$_read_pos;
		}elsif($strand == -1)
		{
			$_dist=$_read_pos-$q_start;
		}
		if($mq==0) { $pList->{'MQ0'}++; }
		if($qscore<$opt_BQ || $mq<$opt_MQ) { next; }
		#if($isTumor)
		#{
		#	if($qscore<$opt_BQ || $mq<$opt_MQ) { next; }
		#}else
		#{
		#	if($qscore<$opt_BQ || $mq<$opt_MQ) { next; }
		#}
		###
		if("\U$qbase" eq "\U$refbase") 
		{ 
			#$_refDepth++; 
			$pList->{'refDepth'}++;
			#push @{$pList->{'refBaseQual'}},$qscore;
			#push @{$pList->{'refMQ'}},$mq;
			#push @{$pList->{'refStrand'}},$strand;
			#push @{$pList->{'base4FET'}},1;
		}else
		{
			$pList->{'ML_ALT_Depth'}++;
			#push @_base,$qbase;
			#push @_bq,$qscore;
			#push @_mq,$mq;
			#push @_d,$_dist;
			#push @_strand,$strand;
		}
		if($out) { printf $out "\tID:$readID\tMQ:$mq\tBQ:$qscore\tread pos:$_read_pos\tbase:$qbase\tstrand:$strand\tdistToEnd:$_dist\tcigar:$cigar\n"; }
	}
}
sub usage
{
	die `pod2text $0`;
}
