#!/usr/bin/env perl
#============================================================================
# Name        		: pgenome_rna.pl 
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Wed Feb 27 11:58:03 2013
# Last Modified By	: 
# Last Modified On	: Wed Feb 27 11:58:03 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl pgenome_rna.pl [option] <config file> <outfile>

	-o		result dir [default: ./analysis]
	-t		number of threads [bwa parameters, default 4]
	-queue		queue [default: ""]
	-project	project [default: ""]
	-ini		ini file [default: $Bin/../parameter/init_human.sh]
	-h		display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime ceil);
use Cwd 'abs_path';
use File::Path;
use File::Basename;
use Graph;
use Graph::Traversal::BFS;
use Data::Dumper;
use FindBin qw($Bin);
#$Data::Dumper::Indent=1;
#$Data::Dumper::Sortkeys=1;

my ($in,$out,$gOut);
my ($opt_h,$opt_o,$opt_queue,$opt_project,$opt_ini,$opt_t);

GetOptions("h"	=>\$opt_h,
	"o=s"=>\$opt_o,
	"queue=s"=>\$opt_queue,
	"project=s"=>\$opt_project,
	"ini=s"=>\$opt_ini,
	"t=i"=>\$opt_t
	);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $outfile=shift @ARGV;
my $g_infile=$infile;
if(!defined($opt_o)) { $opt_o="./analysis"; }
if(!defined($opt_t)) { $opt_t=4; }
if(!defined($opt_queue)) { $opt_queue=""; }
if(!defined($opt_project)) { $opt_project=""; }
if(!defined($opt_ini)) { $opt_ini="$Bin/../parameter/init_human.sh"; } 
$opt_ini="-c $opt_ini";

$opt_o=abs_path($opt_o);

my %config=();
readConfig(\%config,$infile);
open $gOut,">",$outfile or die "$!";

my $g = Graph->new(refvertexed=>1);

my @g_ToDelFile=();

my $pJob={
		name=>"report_One",
		memory=>"500M",
		cmd=>["echo All done!"],
		out=>["$opt_o/report/report.pdf"],
		in=>{},
};

## pipeline specification
foreach my $patientID (keys %config)
{

	#### alignment
	my ($ID_N,$pN);
	$ID_N=(keys %{$config{$patientID}->{'N'}})[0];
	$pN=$config{$patientID}->{'N'}->{$ID_N};

	my $aJob=aln($g,$pN,"$opt_o/$patientID/$ID_N",$ID_N);
	$g->add_edge($aJob,$pJob);

}

genReport($pJob);

## output graph
my $b=Graph::Traversal::BFS->new($g,pre=>\&outputNode);
$b->bfs($g);
my @e=$g->edges;
foreach my $_e (@e) 
{
	outputEdge($_e->[0],$_e->[1]);
}
mkpath("$opt_o/log");
printf $gOut "log_dir $opt_o/log\n";

## delete tmep files
open OUT,">","$opt_o/delete.sh" or die "$!";
foreach (@g_ToDelFile)
{
	printf OUT "rm $_\n";
}
close OUT;

############################################################################


sub aln
{
	my ($g,$pList,$outDir,$ID)=@_;
	my @j=();
	my $cmdFq1="";
	my $cmdFq2="";
	my $cmdFq2_do=1;
	my @mFq1=();
	my @mFq2=();
	foreach my $lib (keys %$pList)
	{
		foreach my $lane (keys %{$pList->{$lib}})
		{
			my $fq1=$pList->{$lib}->{$lane}->[0];
			my $fq2=$pList->{$lib}->{$lane}->[1];
			
			my $job;
			if($fq1 =~ /\.(fq|fastaq|fastq)(\.gz)*$/)
			{
				push @mFq1,$fq1;
			}
			if($fq2 =~ /\.(fq|fastaq|fastq)(\.gz)*$/ && $cmdFq2_do)
			{
				push @mFq2,$fq2;
			}elsif($fq2 eq "-")
			{
				$cmdFq2="echo Concatenation of fq of $ID done.";
				$cmdFq2_do=0;
			}
		}
	}
	my $mergedFq1="$outDir/$ID.R1.fq.gz";
	my $mergedFq2="-";
	if(@mFq1 > 1) { $cmdFq1="cat ". join(" ",@mFq1). " > $mergedFq1"; }
	else { $cmdFq1="ln -s $mFq1[0] $mergedFq1"; }
	push @g_ToDelFile,$mergedFq1;
	if($cmdFq2_do)
	{
		$mergedFq2="$outDir/$ID.R2.fq.gz";
		if(@mFq2 > 1) { $cmdFq2="cat ". join(" ",@mFq2). " > $mergedFq2"; }
		else { $cmdFq2="ln -s $mFq2[0] $mergedFq2"; }
		push @g_ToDelFile,$mergedFq2;
	}
	#merge
	my $mJob={
		name=>"merge_fastaq_$ID",
		memory=>"1G",
		cmd=>[$cmdFq1,$cmdFq2],
		out=>[$mergedFq1,$mergedFq2],
	};
	$g->add_vertex($mJob);

	#aln & count reads & variants
	my $aJob={
		name=>"aln_count_var_$ID",
		memory=>"80G",
		queue=>"mem.q",
		project=>"mem",
		cpuNum=>16,
		cmd=>["$Bin/aln.hisat.sh -t 12 $outDir/hisat $mergedFq1 $mergedFq2 $ID unknown unknown"],
		out=>["$outDir/hisat/$ID.hisat.hit.sort.bam"],
	};
	$g->add_edge($mJob,$aJob);
	
	return ($aJob);
}

sub genReport
{
	my ($pJob)=@_;

}

sub readConfig
{
	my ($pList,$infile)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		if(@F<9) { usage(); }
		#if(!exists($pList->{$F[0]})) { $pList->{$F[0]}={"N"=>{'sampleID'=>'','data'=>{}},"T"=>{'sampleID'=>'','data'=>{}}}; }

		#$pList->{$F[0]}->{$F[1]}->{'sampleID'}=$F[2];
		#$pList->{$F[0]}->{$F[1]}->{'data'}->{"$F[3]"}->{"$F[4]"}=[$F[5],$F[6]];
		
		####
		$F[0]=~s/ /_/g;
		$F[0]=~s/:/-/g;
		$F[0]=~s/,/-/g;
		$F[0]=~s/"//g;
		$F[0]=~s/\//-/g;
		$F[2]=~s/ /_/g;
		$F[2]=~s/:/-/g;
		$F[2]=~s/,/-/g;
		$F[2]=~s/"//g;
		$F[2]=~s/\//-/g;
		if(!exists($pList->{$F[0]})) { $pList->{$F[0]}={"N"=>{},"T"=>{}}; }
		$pList->{$F[0]}->{$F[1]}->{$F[2]}->{$F[3]}->{$F[4]}=[$F[5],$F[6],$F[7],$F[8]];

		#patientID   sampleType  sampleID    libID   lane    fq1 fq2	fragmentCount	qcStat	other
		#A549    N       A549    LIB2501 L7      /Bak/admin/cell_line675/fq/R141/LIB2501_SAM635619_L7_R1.fastq.gz        /Bak/admin/cell_line675/fq/R141/LIB2501_SAM635619_L7_R2.fastq.gz        0       0       Lung
		
	}
}



sub outputEdge
{
	my ($u,$v,$b)=@_;
	printf $gOut "order %s before %s\n",$u->{'name'},$v->{'name'};
	
}
sub outputNode
{
	my ($v,$b)=@_;
	#printf STDERR "%s\t%s\n",$v,$v->{'name'};
	my $_cmd=$v->{'cmd'};
	my $_name=$v->{'name'};
	
	## output shell
	#printf STDERR "$v->{'out'}->[0]\n";
	my $_outDir=dirname($v->{'out'}->[0]);
	mkpath($_outDir);
	open $out,">","$_outDir/$_name.sh" or die "Die in generate $_outDir/$_name.sh ($!)";
	printf $out "%s\n",join(" && \\\n",@$_cmd);
	close $out;

	## output sjm job
	printf $gOut "job_begin\n";
	printf $gOut "  name $v->{name}\n";
	printf $gOut "  memory $v->{memory}\n";
	if(defined($v->{"cpuNum"})) { printf $gOut "  cpu $v->{cpuNum}\n"; }
	if(defined($v->{"queue"})) { printf $gOut "  queue $v->{queue}\n"; }
	if(defined($v->{"project"})) { printf $gOut "  project $v->{project}\n"; }
	#if($opt_queue ne "") { printf $gOut "  queue $opt_queue\n"; }
	#if($opt_project ne "") { printf $gOut "  project $opt_project\n"; }
	printf $gOut "  cmd_begin\n";
	@$_cmd=map { "    $_" } @$_cmd;
	printf $gOut "%s\n",join(" &&\n",@$_cmd);
	printf $gOut "  cmd_end\n";
	printf $gOut "job_end\n";
}
sub usage
{
	die `pod2text $0`;
}
