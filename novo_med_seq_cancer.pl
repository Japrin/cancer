#!/usr/bin/env perl
#============================================================================
# Name        		: novo_wes_cancer.pl
# Author      		: zhengliangtao (tao2013@gmail.com)
# Version     		: v1.00
# Created On  		: Wed Feb 27 11:58:03 2013
# Last Modified By	: 
# Last Modified On	: Wed Feb 27 11:58:03 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl novo_med_seq_cancer.pl [option] <config file> <outfile>

	-o		result dir [default: ./analysis]
	-p		extra 'aln.mem.sh' parameter [quoted by "", default ""]
	-t		number of threads [bwa parameters, default 4]
	-javaM		java(picard, gatk etc) memory requirement, ex. 32G [default ""]
	-alnM		bwa memory requirement, ex. 32G [default ""]
	-normdup	don't remove duplicate [default: OFF]
	-realn		do re-alignemnt [default: OFF]
	-recal		do base quality recalibaration [default: OFF]
	-var		call snp/indel per sample [default: OFF]
	-single		single sample (must be 'N' sampleType) [default: OFF]
	-nosomatic	do call somatic variation [default: OFF]
	-TR		Target Region file [default: ""]
	-binSize	binSize for single sample WGS CNV calling [default 100]
	-cna		do somatic CNV analysis [default: OFF]
	-ssv		do somatic SV analysis [default: OFF]
	-queue		queue [default: ""]
	-project	project [default: ""]
	-maxreads	max reads in fq file; if specified, must >= 1000,000]
	-bychr		post-alignment process, snp/indel calling etc by chromosomes in file 'bychr' [default OFF ]
	-ini		ini file [default: $Bin/parameter/init_human.sh]
	-samtools	using samtools to do snp&indel calling [default OFF]
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
my ($opt_h,$opt_bin,$opt_o,$opt_m,$opt_realn,$opt_recal,$opt_normdup,$opt_nosomatic,$opt_var,$opt_TR,$opt_queue,$opt_project,$opt_ini,$opt_cna,$opt_single,$opt_ssv,$opt_p,$opt_maxreads,$opt_bychr,$opt_t,$opt_javaM,$opt_alnM,$opt_binSize,$opt_samtools);
GetOptions("h"	=>\$opt_h,"bin=s"=>\$opt_bin,"o=s"=>\$opt_o,"realn"=>\$opt_realn,"recal"=>\$opt_recal,"normdup"=>\$opt_normdup,"nosomatic"=>\$opt_nosomatic,"var"=>\$opt_var,"TR=s"=>\$opt_TR,"queue=s"=>\$opt_queue,"project=s"=>\$opt_project,"ini=s"=>\$opt_ini,"cna"=>\$opt_cna,"single"=>\$opt_single,"ssv"=>\$opt_ssv,"p=s"=>\$opt_p,"maxreads=i"=>\$opt_maxreads,"bychr=s"=>\$opt_bychr,"t=i"=>\$opt_t,"m=s"=>\$opt_m,"javaM=s"=>\$opt_javaM,"alnM=s"=>\$opt_alnM,"binSize=i"=>\$opt_binSize,"samtools"=>\$opt_samtools);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $outfile=shift @ARGV;
my $g_infile=$infile;
if(!defined($opt_bin)) { $opt_bin=$Bin; }
if(!defined($opt_o)) { $opt_o="./analysis"; }
if(!defined($opt_queue)) { $opt_queue=""; }
if(!defined($opt_project)) { $opt_project=""; }
if(!defined($opt_TR)) { $opt_TR=""; } else { $opt_TR="-r $opt_TR"; }
if(!defined($opt_ini)) { $opt_ini="$Bin/parameter/init_human.sh"; } 
$opt_ini="-c $opt_ini";
if(!defined($opt_p)) { $opt_p=""; } else { $opt_p="-p \"$opt_p\""; }
if(defined($opt_maxreads)) { $opt_maxreads=$opt_maxreads >= 1000000 ? $opt_maxreads : 1000000; }
if(!defined($opt_t)) { $opt_t=4; }
if(!defined($opt_m)) { $opt_m=""; }
if(!defined($opt_javaM)) { $opt_javaM=""; }
if(!defined($opt_alnM)) { $opt_alnM=""; }
if(!defined($opt_binSize)) { $opt_binSize=100; }
#if(!defined($opt_bychr)) { $opt_bychr="/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/chr.list"; }

my $gmm=($opt_javaM eq ""?"16G":"$opt_javaM");
$gmm=~s/g//i;
$opt_o=abs_path($opt_o);

my %config=();
readConfig(\%config,$infile);
open $gOut,">",$outfile or die "$!";

my $g = Graph->new(refvertexed=>1);

my @g_ToDelFile=();

my $pJob={
		name=>"report_One",
		memory=>"500M",
		cmd=>[],
		out=>["$opt_o/report/report.pdf"],
		in=>{},
};
foreach my $patientID (keys %config)
{

	#### alignment
	my ($ID_N,$pN);
	$ID_N=(keys %{$config{$patientID}->{'N'}})[0];
	$pN=$config{$patientID}->{'N'}->{$ID_N};
	my ($bamN,$bamNOne,$jAlnN,$jFinalN,$chrN)=aln($g,$pN,"$opt_o/$patientID/$ID_N/aln",$ID_N);
	my @jArrayVarN=();
	my @jArrayVarT=();
	
	my @alnJob=();  ## element is array reference: [$bamT,$bamTOne,jAlnT,chrT]

	if(!defined($opt_single))
	{
		my @TSamples=keys %{$config{$patientID}->{'T'}};
		foreach my $ID_T (@TSamples)
		{
			my $pT=$config{$patientID}->{'T'}->{$ID_T};
			my ($bamT,$bamTOne,$jAlnT,$jFinalT,$chrT)=aln($g,$pT,"$opt_o/$patientID/$ID_T/aln",$ID_T);
			push @alnJob,[$ID_T,$bamT,$bamTOne,$jAlnT,$jFinalT,$chrT];
		}
	}
	
	#### variation in every sample
	if($opt_var)
	{
		### 'N' sample
		my %_array_N=();
		for(my $i=0;$i<@$bamN;$i++)
		{
			my %jVarN=var($g,$bamN->[$i],"$opt_o/$patientID/$ID_N/var/byChr/$chrN->[$i]","$ID_N.$chrN->[$i]");
			foreach (keys %jVarN)
			{
				$g->add_edge($jAlnN->[$i],$jVarN{$_});
				push @{$_array_N{$_}},$jVarN{$_};
			}
		}
		foreach (keys %_array_N)
		{
			my $mergeVarJob=mergeVar($g,$_array_N{$_},"$patientID.$ID_N",$_,"$opt_o/$patientID/$ID_N/var");
			$g->add_edge($mergeVarJob,$pJob);
			$pJob->{'in'}->{'var'}->{$_}->{$ID_N}=$mergeVarJob;
		}
		### variation by whole genome
		my %jVar=varOne($g,$bamNOne,"$opt_o/$patientID/$ID_N/sv","$patientID.$ID_N");
		foreach (keys %jVar)
		{
			$g->add_edge($jFinalN,$jVar{$_});
			$g->add_edge($jVar{$_},$pJob);
			$pJob->{'in'}->{'var'}->{$_}->{$ID_N}=$jVar{$_};
		}
		### 'T' sample
		if(!defined($opt_single))
		{
			foreach my $_s (@alnJob)
			{
				my ($ID_T,$bamT,$bamTOne,$jAlnT,$jFinalT,$chrT)=@$_s;
				### variation by chr
				my %_array_T=();
				for(my $i=0;$i<@$bamT;$i++)
				{
					my %jVarT=var($g,$bamT->[$i],"$opt_o/$patientID/$ID_T/var/byChr/$chrT->[$i]","$ID_T.$chrT->[$i]");
					foreach (keys %jVarT)
					{
						$g->add_edge($jAlnT->[$i],$jVarT{$_});
						push @{$_array_T{$_}},$jVarT{$_};
					}
				}
				foreach (keys %_array_T)
				{
					my $mergeVarJob=mergeVar($g,$_array_T{$_},"$patientID.$ID_T",$_,"$opt_o/$patientID/$ID_T/var");
					$g->add_edge($mergeVarJob,$pJob);
					$pJob->{'in'}->{'var'}->{$_}->{$ID_T}=$mergeVarJob;
				}
				### variation by whole genome
				my %jVar=varOne($g,$bamTOne,"$opt_o/$patientID/$ID_T/sv","$patientID.$ID_T");
				foreach (keys %jVar)
				{
					$g->add_edge($jFinalT,$jVar{$_});
					$g->add_edge($jVar{$_},$pJob);
					$pJob->{'in'}->{'var'}->{$_}->{$ID_T}=$jVar{$_};
				}
			}
		}
	}
	#### coverage
	my $jStatN=statAln($g,$bamNOne,"$opt_o/$patientID/$ID_N/aln/stat","$ID_N");
	$g->add_edge($jFinalN,$jStatN);
	$g->add_edge($jStatN,$pJob);
	$pJob->{'in'}->{'stat'}->{$ID_N}=$jStatN;
	#my @_array_N=();
	#for(my $i=0;$i<@$bamN;$i++)
	#{
	#	my $jStatN=statAln($g,$bamN->[$i],"$opt_o/$patientID/$ID_N/aln/stat/byChr/$chrN->[$i]","$ID_N.$chrN->[$i]");
	#	$g->add_edge($jAlnN->[$i],$jStatN);
	#	push @_array_N,$jStatN;
	#}
	#my $mergeStatJob=mergeStat($g,\@_array_N,"$patientID.$ID_N","$opt_o/$patientID/$ID_N/aln/stat");
	#$g->add_edge($mergeStatJob,$pJob);
	#$pJob->{'in'}->{'stat'}->{$ID_N}=$mergeStatJob;
	if(!defined($opt_single))
	{
		foreach my $_s (@alnJob)
		{
			my ($ID_T,$bamT,$bamTOne,$jAlnT,$jFinalT,$chrT)=@$_s;
			my $jStatT=statAln($g,$bamTOne,"$opt_o/$patientID/$ID_T/aln/stat","$ID_T");
			$g->add_edge($jFinalT,$jStatT);
			$g->add_edge($jStatT,$pJob);
			$pJob->{'in'}->{'stat'}->{$ID_T}=$jStatT;


			#my @_array_T=();
			#for(my $i=0;$i<@$bamT;$i++)
			#{
			#	my $jStatT=statAln($g,$bamT->[$i],"$opt_o/$patientID/$ID_T/aln/stat/byChr/$chrT->[$i]","$ID_T.$chrT->[$i]");
			#	$g->add_edge($jAlnT->[$i],$jStatT);
			#	push @_array_T,$jStatT;
			#}
			#my $mergeStatJob=mergeStat($g,\@_array_T,"$patientID.$ID_T","$opt_o/$patientID/$ID_T/aln/stat");
			#$g->add_edge($mergeStatJob,$pJob);
			#$pJob->{'in'}->{'stat'}->{$ID_T}=$mergeStatJob;
		}
	}
	
	#### somatic
	if(!$opt_nosomatic && !defined($opt_single))
	{
		foreach my $_s (@alnJob)
		{
			my ($ID_T,$bamT,$bamTOne,$jAlnT,$jFinalT,$chrT)=@$_s;
			### calling by chr
			my %_array=();
			for(my $i=0;$i<@$bamT;$i++)
			{
				my %jSomatic=somatic($g,$bamN->[$i],$bamT->[$i],"$opt_o/$patientID/somatic/$ID_T/byChr/$chrT->[$i]","$patientID.$ID_T.$chrT->[$i]");
				foreach (keys %jSomatic)
				{
					$g->add_edge($jAlnN->[$i],$jSomatic{$_});
					$g->add_edge($jAlnT->[$i],$jSomatic{$_});
					push @{$_array{$_}},$jSomatic{$_};
				}
			}
			foreach (keys %_array)
			{
				my $mergeVarJob=mergeVar($g,$_array{$_},"$patientID.$ID_T","$_","$opt_o/$patientID/somatic/$ID_T/$_","somatic","$patientID.$ID_N");
				$g->add_edge($mergeVarJob,$pJob);
				$pJob->{'in'}->{'somatic'}->{$_}->{$ID_T}=$mergeVarJob;
			}
			### calling as a whole
			my %jSomatic=somaticOne($g,$bamNOne,$bamTOne,"$opt_o/$patientID/somatic/$ID_T","$patientID.$ID_T",$pJob->{'in'}->{'var'}->{"sample_GATK"}->{$ID_N});
			foreach (keys %jSomatic)
			{
				$g->add_edge($jFinalT,$jSomatic{$_});
				$g->add_edge($jFinalN,$jSomatic{$_});
				$g->add_edge($jSomatic{$_},$pJob);
				$pJob->{'in'}->{'somatic'}->{$_}->{$ID_T}=$jSomatic{$_};
			}
		}
	}
}
genReport($pJob);

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
		if(@F<9) { print STDERR "config file format error!\n"; usage(); }
		#if(!exists($pList->{$F[0]})) { $pList->{$F[0]}={"N"=>{'sampleID'=>'','data'=>{}},"T"=>{'sampleID'=>'','data'=>{}}}; }

		#$pList->{$F[0]}->{$F[1]}->{'sampleID'}=$F[2];
		#$pList->{$F[0]}->{$F[1]}->{'data'}->{"$F[3]"}->{"$F[4]"}=[$F[5],$F[6]];
		
		####
		if(!exists($pList->{$F[0]})) { $pList->{$F[0]}={"N"=>{},"T"=>{}}; }
		$pList->{$F[0]}->{$F[1]}->{$F[2]}->{$F[3]}->{$F[4]}=[$F[5],$F[6],$F[7],$F[8]];

		#patientID   sampleType  sampleID    libID   lane    fq1 fq2	fragmentCount	qcStat
		
	}
}
sub aln
{
	my ($g,$pList,$outDir,$ID)=@_;
	my @j=();
	foreach my $lib (keys %$pList)
	{
		foreach my $lane (keys %{$pList->{$lib}})
		{
			my $fq1=$pList->{$lib}->{$lane}->[0];
			my $fq2=$pList->{$lib}->{$lane}->[1];
			my $nRead=$pList->{$lib}->{$lane}->[2];
			
			my $job;
			if($fq1 =~ /\.(fq|fastaq|fastq)(\.gz)*$/)
			{
				if(defined($opt_maxreads) && $nRead>$opt_maxreads)
				{
					my $_n=ceil($nRead/$opt_maxreads);
					my $_prefix1="$outDir/".basename($fq1);
					my $_prefix2="$outDir/".basename($fq2);
					$_prefix1 =~ s/\.fq(\.gz)*$//;
					$_prefix2 =~ s/\.fq(\.gz)*$//;

					### split fq already done
					my $_fff=1;
					for(my $i=1;$i<=$_n;$i++)
					{
						my $_SFq1=sprintf "$_prefix1.S%06d.fq.gz",$i;
						my $_SFq2=sprintf "$_prefix2.S%06d.fq.gz",$i;
						if(!(-e $_SFq1) || !(-e $_SFq2))
						{
							$_fff=0;
							last;
						}
					}
					my $_ccc1="";
					my $_ccc2="";
					if($_fff)
					{
						#don't split fq
						$_ccc1="echo \"already done: $fq1\"";
						$_ccc2="echo \"already done: $fq2\"";
					}else
					{
						#split fq
						$_ccc1="$opt_bin/alignment/split.fq.pl -n $opt_maxreads $fq1 $_prefix1";
						$_ccc2="$opt_bin/alignment/split.fq.pl -n $opt_maxreads $fq2 $_prefix2";
					}

					my $_job0={
						name=>"split_fq_${lib}_${lane}",
						memory=>"500m",
						cmd=>[$_ccc1,$_ccc2],
						#cmd=>["$opt_bin/alignment/split.fq.pl -n $opt_maxreads $fq1 $_prefix1","$opt_bin/alignment/split.fq.pl -n $opt_maxreads $fq2 $_prefix2"],
						out=>["$_prefix1","$_prefix2"],
					};
					$g->add_vertex($_job0);
					my @s=();
					for(my $i=1;$i<=$_n;$i++)
					{
						my $_mmOpt=($opt_javaM eq ""?"":"-m $gmm");
						my $_SFq1=sprintf "$_prefix1.S%06d.fq.gz",$i;
						my $_SFq2=sprintf "$_prefix2.S%06d.fq.gz",$i;
						my $_job={
							name=>"aln_${lib}_${lane}_${i}",
							cpuNum=>$opt_t,
							memory=>($opt_alnM eq ""?"16G":"$opt_alnM"),
							cmd=>["$opt_bin/alignment/aln.mem.sh -i $i -t $opt_t $_mmOpt $opt_p $opt_ini $outDir $_SFq1 $_SFq2 $ID $lib $lane"],
							out=>["$outDir/$ID.$lib.$lane.$i.sort.bam"],
						};
						$g->add_edge($_job0,$_job);
						push @s,$_job;
						### to be deleted
						push @g_ToDelFile,$_SFq1,$_SFq2,"$outDir/$ID.$lib.$lane.$i.sort.bam";
					}
	
					my $_bams="";
					foreach (@s){ $_bams.=" $_->{'out'}->[0]"; }
					my $_outBam="$outDir/$ID.$lib.$lane.sort.bam";
	
					$job={
						name=>"primary_aln_${ID}_${lib}_${lane}",
						memory=>(@s>1)?"7500M":"100M",
						cmd=>["$opt_bin/alignment/merge_bam.sh $opt_ini $_outBam $_bams"],
						out=>[$_outBam],
					};
					foreach (@s) { $g->add_edge($_,$job); }
					### to be deleted
					push @g_ToDelFile,$_outBam;

				}else
				{
					my $_mmOpt=($opt_javaM eq ""?"":"-m $gmm");
					$job={
						name=>"primary_aln_${ID}_${lib}_${lane}",
						memory=>($opt_alnM eq ""?"16G":"$opt_alnM"),
						cpuNum=>$opt_t,
						cmd=>["$opt_bin/alignment/aln.mem.sh -t $opt_t $_mmOpt $opt_p $opt_ini $outDir $fq1 $fq2 $ID $lib $lane"],
						out=>["$outDir/$ID.$lib.$lane.sort.bam"],
					};
					$g->add_vertex($job);
					### to be deleted
					push @g_ToDelFile,"$outDir/$ID.$lib.$lane.sort.bam";
				}
			}else
			{
				$job={
					name=>"primary_aln_${ID}_${lib}_${lane}",
					memory=>"100m",
					cmd=>["echo \"already done: $fq1\""],
					out=>["$fq1"],
				};
				$g->add_vertex($job);
			}
			push @j,$job;
		}
	}
	#merge
	my $bams="";
	foreach (@j){ $bams.=" $_->{'out'}->[0]"; }
	my $retBam="$outDir/$ID.bam";	
	my $_mmOpt=($opt_javaM eq ""?"":"-m $gmm");
	my $mJob={
		name=>"merge_bam_primary_$ID",
		memory=>(@j>1)?($opt_javaM eq ""?"10000M":$opt_javaM):"100M",
		cmd=>["$opt_bin/alignment/merge_bam.sh $opt_ini $_mmOpt $retBam $bams"],
		out=>["$retBam"],
	};
	### to be deleted
	#push @g_ToDelFile,"$retBam";
	
	if(!$opt_normdup)
	{
		my $_bam=$retBam;
		$retBam=~s/\.bam$/\.nodup\.bam/;
		push @{$mJob->{'cmd'}},"$opt_bin/alignment/rmdup.sh $opt_ini $_mmOpt $_bam $retBam";
		$mJob->{'memory'}=($opt_javaM eq ""?"10000M":$opt_javaM);
		$mJob->{'out'}->[0]=$retBam;
		$mJob->{'name'}="merge_markDup_bam_all_$ID";
		### to be deleted
		push @g_ToDelFile,"$_bam";
	}
	$g->add_vertex($mJob);
	foreach (@j) { $g->add_edge($_,$mJob); }
	### byChr
	my @fJob=();
	my @fChr=();
	if($opt_bychr && (-e $opt_bychr))
	{
		open IN,$opt_bychr or die "$!";
		my @chrList=<IN>;
		chomp @chrList;
		close IN;
		foreach my $_chr (@chrList)
		{
			my $_ifile=$mJob->{'out'}->[0];
			my $_bname=basename($_ifile);
			$_bname=~s/\.bam$/\.$_chr\.bam/;
			my $_ofile=dirname($_ifile)."/byChr/$_chr/".$_bname;
			my $_chrJob={
				name=>"extract_bam_${_chr}_$ID",
				memory=>"4500M",
				cmd=>["$opt_bin/alignment/bin_sam.sh $_chr $_ofile $_ifile"],
				out=>["$_ofile"],
			};
			$g->add_edge($mJob,$_chrJob);
			push @fJob,$_chrJob;
			push @fChr,$_chr;

		}
		### to be deleted
		push @g_ToDelFile,$mJob->{'out'}->[0];
	}else
	{
		push @fJob,$mJob;
		push @fChr,"all";
	}
	
	### realn
	if($opt_realn)
	{
		foreach my $_chrJob (@fJob)
		{
			my $_iBam=$_chrJob->{'out'}->[0];
			my $_oBam=$_iBam;
			$_oBam=~s/\.bam$/\.realn\.bam/;
			push @{$_chrJob->{'cmd'}},"$opt_bin/alignment/realn.sh $opt_ini $_iBam $_oBam";
			$_chrJob->{'memory'}="10000M";
			$_chrJob->{'out'}->[0]=$_oBam;
			my ($_title,$_chr)=$_chrJob->{'name'}=~/(.+?)_bam_(.+)/;
			$_chrJob->{'name'}="${_title}_realn_bam_${_chr}";
			### to be deleted
			push @g_ToDelFile,$_iBam;
		}
	}
	### recalibration
	if($opt_recal)
	{
		foreach my $_chrJob (@fJob)
		{
			my $_iBam=$_chrJob->{'out'}->[0];
			my $_oBam=$_iBam;
			$_oBam=~s/\.bam$/\.recal\.bam/;
			push @{$_chrJob->{'cmd'}},"$opt_bin/alignment/recal.sh $opt_ini $_iBam $_oBam";
			$_chrJob->{'memory'}="10000M";
			$_chrJob->{'out'}->[0]=$_oBam;
			my ($_title,$_chr)=$_chrJob->{'name'}=~/(.+?)_bam_(.+)/;
			$_chrJob->{'name'}="${_title}_recal_bam_${_chr}";
			### to be deleted
			push @g_ToDelFile,$_iBam;
		}
	}
	my @retBam=();
	foreach my $_job (@fJob) { push @retBam,$_job->{'out'}->[0]; }
	
	### final merged one bam
	my $finalBam="$outDir/$ID.final.bam";
	$bams=join(" ",@retBam);
	my $finalMergeJob={
		name=>"merge_bam_final_$ID",
		memory=>"7500M",
		cmd=>["$opt_bin/alignment/merge_bam.sh $opt_ini $finalBam $bams","md5sum $finalBam > $finalBam.md5"],
		out=>["$finalBam"],
	};
	foreach my $_job (@fJob) { $g->add_edge($_job,$finalMergeJob); }

	return (\@retBam,$finalBam,\@fJob,$finalMergeJob,\@fChr);
}
sub statAln
{
	my ($g,$bam,$outDir,$ID)=@_;
	my $mJob={
		name=>"stat_$ID",
		memory=>"5G",
		cmd=>["$opt_bin/stat/stat.aln.sh $opt_ini $opt_TR $bam $outDir $ID"],
		out=>["$outDir/${ID}_mapping_coverage.txt","$outDir/depth_frequency.xls","$outDir/cumu.xls","$outDir/$ID.coverage.bychr.txt"],
	};
	return $mJob;
}
sub mergeStat
{
	my ($g,$pJobArray,$ID,$outDir)=@_;
	my $_job={
		name=>"merge_stat_${ID}",
		memory=>"500m",
		cmd=>[],
		out=>[],
	};
	## table
	my $infile="";
	foreach (@$pJobArray)
	{
		$infile.=" $_->{'out'}->[0]";
		$g->add_edge($_,$_job);
	}
	my $afile="$outDir/$ID.coverge.stat.final.txt";
	push @{$_job->{'cmd'}},"$opt_bin/stat/stat.merge.sh $afile $infile";
	push @{$_job->{'out'}},$afile;
	## depth
	$infile="";
	foreach (@$pJobArray) { $infile.=" $_->{'out'}->[1]"; }
	$afile="$outDir/$ID.depth.txt";
	push @{$_job->{'cmd'}},"$opt_bin/stat/stat.merge.sh $afile $infile";
	push @{$_job->{'out'}},$afile;
	## cumu depth
	$infile="";
	foreach (@$pJobArray) { $infile.=" $_->{'out'}->[2]"; }
	$afile="$outDir/$ID.cumu.txt";
	push @{$_job->{'cmd'}},"$opt_bin/stat/stat.merge.sh $afile $infile";
	push @{$_job->{'out'}},$afile;
	## coverage by chr
	$infile="";
	foreach (@$pJobArray) { $infile.=" $_->{'out'}->[3]"; }
	$afile="$outDir/$ID.covByChr.txt";
	push @{$_job->{'cmd'}},"$opt_bin/stat/stat.merge.sh $afile $infile";
	push @{$_job->{'out'}},$afile;


	return $_job;
}

sub var
{
	my ($g,$bam,$outDir,$ID)=@_;
	my %ret=();
	my $mJob;
	if($opt_TR)
	{
		if($opt_samtools)
		{
			$mJob={
				name=>"var_$ID",
				memory=>"2G",
				cmd=>["$opt_bin/var/var_samtools.sh $opt_ini $opt_TR $ID $bam $outDir"],
				out=>["$outDir/$ID.samtools.snp.vcf","$outDir/$ID.samtools.indel.vcf"],
				type=>["snp","indel"],
			};
		}else
		{
			$mJob={
				name=>"var_$ID",
				memory=>"7G",
				cmd=>["$opt_bin/var/var_gatk.sh $opt_ini $opt_TR $ID $bam $outDir"],
				out=>["$outDir/$ID.GATK.snp.flt.vcf","$outDir/$ID.GATK.indel.flt.vcf"],
				type=>["snp","indel"],
			};
		}
	}else
	{
		if($opt_samtools)
		{
			$mJob={
				name=>"var_$ID",
				memory=>"7G",
				cmd=>["$opt_bin/var/var_samtools.sh $opt_ini $ID $bam $outDir"],
				out=>["$outDir/$ID.samtools.snp.vcf","$outDir/$ID.samtools.indel.vcf"],
				type=>["snp","indel"],
			};
		}else
		{
			$mJob={
				name=>"var_$ID",
				memory=>"7G",
				cmd=>["$opt_bin/var/var_gatk.sh $opt_ini $ID $bam $outDir"],
				out=>["$outDir/$ID.GATK.snp.flt.vcf","$outDir/$ID.GATK.indel.flt.vcf"],
				type=>["snp","indel"],
				#cmd=>["$opt_bin/var/var_gatk_haplotypeCaller.sh $opt_ini $ID $bam $outDir"],
				#out=>["$outDir/$ID.GATK.var.raw.vcf"],
				#type=>["snp_indel"],
			};
		}
	}
	$ret{'sample_GATK'}=$mJob;
	return %ret;
}
sub varOne
{
	my ($g,$bam,$outDir,$ID)=@_;
	my %ret=();
	if(!$opt_TR)
	{
		my $_job={
			name=>"sample_breakdancer_$ID",
			memory=>"7G",
			cmd=>["$opt_bin/var/var_sv_breakdancer.sh $opt_ini $ID $outDir/breakdancer $bam"],
			out=>["$outDir/breakdancer/$ID.breakdancer.gff"],
		};
		$ret{'sample_breakdancer'}=$_job;

		#my $_job1={
		#	name=>"sample_pindel_$ID",
		#	memory=>"7G",
		#	cmd=>["$opt_bin/var/var_sv_pindel.sh -b $outDir/breakdancer/$ID.breakdancer.txt $opt_ini $ID $outDir/pindel $bam"],
		#	out=>["$outDir/pindel/$ID.pindel.gff"],
		#};
		#$ret{'sample_pindel'}=$_job1;
		#$g->add_edge($_job,$_job1);

		#cnv only for whole genome
		my $_job2={
			name=>"sample_cnvnator_$ID",
			memory=>"7G",
			cmd=>["$opt_bin/var/var_sv_cnvnator.sh $opt_ini -b $opt_binSize $ID $outDir/cnvnator $bam"],
			out=>["$outDir/cnvnator/$ID.cnvnator.gff"],
		};
		$ret{'sample_cnvnator'}=$_job2;
	}
	return %ret;
}
sub mergeVar
{
	my ($g,$pJobArray,$ID,$method,$outDir,$pa,$ID2)=@_;
	if(!defined($pa)) { $pa="\"\""; }
	my $aMethod="";
	if($method =~ /GATK/)
	{
		$aMethod="GATK";
	}
	elsif($method ne "GATK" && $method ne "varScan" && $method ne "samtools")
	{
		$aMethod="with.format";
	}else
	{
		$aMethod=$method;
	}
	my $annID;
	if(!defined($ID2)) { $annID="$ID"; }
	else 
	{ 
		#$annID="$ID2,$ID"; 
		$annID="N,T"; 
	}

	my @type=@{$pJobArray->[0]->{'type'}};
	my $_n=@{$pJobArray->[0]->{'out'}};
	my $_m=@$pJobArray;
	my $_job={
		name=>"merge_var_${method}_${ID}",
		memory=>"5G",
		cmd=>[],
		out=>[],
	};
	for(my $i=0;$i<$_n;$i++)
	{
		my $invcf="";
		foreach (@$pJobArray) { $invcf.=" $_->{'out'}->[$i]"; }
		my $avcf="$outDir/$ID.$method.$type[$i].call.vcf";
		push @{$_job->{'cmd'}},"$opt_bin/toolkit/vcf_merge.sh $opt_ini $avcf $invcf";
		if($opt_TR)
		{
			## nothing to do
		}else
		{
			#push @{$_job->{'cmd'}},"$opt_bin/var/var_gatk_varRecalibrate.sh $opt_ini $avcf $avcf.flt.vcf";
			#push @{$_job->{'cmd'}},"mv $avcf.flt.vcf $avcf";
		}

		push @{$_job->{'cmd'}},"$opt_bin/var/var_annotation.sh $opt_ini -m $aMethod -a $pa $avcf $annID";
		push @{$_job->{'out'}},"$outDir/$ID.$method.$type[$i].call.reformated.vcf.gz";
	}
	foreach (@$pJobArray)
	{
		$g->add_edge($_,$_job);
	}
	return $_job;
}
sub somatic
{
	my ($g,$bamN,$bamT,$outDir,$ID)=@_;
	my $mJob_1={
		name=>"somatic_strelka_$ID",
		memory=>"8500M",
		cmd=>["$opt_bin/somatic_snv/somatic_strelka.sh $opt_ini $opt_TR $ID $bamN $bamT $outDir/somatic_strelka"],
		out=>["$outDir/somatic_strelka/results/$ID.strelka.somatic.snvs.filter.vcf","$outDir/somatic_strelka/results/$ID.strelka.somatic.indels.filter.vcf"],
		type=>["snv","indel"],
	};
	my $mJob_2={
		name=>"somatic_varScan_$ID",
		memory=>"6500M",
		cmd=>["$opt_bin/somatic_snv/somatic_varScan.sh $opt_ini $opt_TR $ID $bamN $bamT $outDir/varScan"],
		#out=>["$outDir/varScan/$ID.varScan.call.vcf"],
		out=>["$outDir/varScan/$ID.varScan.call.Somatic.vcf","$outDir/varScan/$ID.varScan.call.Germline.vcf","$outDir/varScan/$ID.varScan.call.LOH.vcf"],
		type=>["somatic.snv","germline.snv","loh.snv"],
	};
	$g->add_vertex($mJob_1);
	$g->add_vertex($mJob_2);
	my %ret=(somatic_strelka=>$mJob_1,varScan=>$mJob_2);
	return %ret;
}
sub somaticOne
{
	my ($g,$bamN,$bamT,$outDir,$ID,$eJob)=@_;
	my %ret=();
	if($opt_cna)
	{
		my $_Job={
			name=>"somatic_cnv_$ID",
			memory=>"6500M",
			cmd=>[],
			out=>[],
			type=>["scnv"],
		};
		if($opt_TR)
		{
			#exome cnv
			#push @{$_Job->{'cmd'}},"$opt_bin/somatic_cna/ecnv_cnv.sh $opt_ini $opt_TR $ID $bamN $bamT $outDir/TR_cnv";
			#push @{$_Job->{'out'}},"$outDir/TR_cnv/$ID.cnv.cnv.txt";
			#push @{$_Job->{'cmd'}},"$opt_bin/somatic_cna/somatic_cnv_varscan.sh $opt_ini $opt_TR $ID $bamN $bamT $outDir/TR_cnv";
			#push @{$_Job->{'out'}},"$outDir/TR_cnv/$ID.varcan.cnv.ann.txt";
			my $normal_snp=$eJob->{'out'}->[0];
			push @{$_Job->{'cmd'}},"$opt_bin/somatic_cna/somatic_cnv_ADTEx.sh -b $normal_snp $opt_ini $opt_TR $ID $bamN $bamT $outDir/TR_cnv";
			push @{$_Job->{'out'}},"$outDir/TR_cnv/out/cnv.result.segment.CNV.gff";
			$g->add_edge($eJob,$_Job);
		}else
		{
			#whole genome cnv
			push @{$_Job->{'cmd'}},"$opt_bin/var/var_sv_freec.sh $opt_ini $ID $outDir/WGS_cnv $bamT $bamN";
			push @{$_Job->{'out'}},"$outDir/WGS_cnv/".basename($bamT)."_CNVs";
		}
		$g->add_vertex($_Job);
		$ret{'somatic_cna'}=$_Job;
	}
	if($opt_ssv)
	{
		my $_Job={
			name=>"somatic_sv_$ID",
			memory=>"6500M",
			cmd=>["$opt_bin/somatic_sv/somatic_breakdancer.sh $opt_ini $ID $bamN $bamT $outDir/breakdancer_ssv"],
			out=>["$outDir/breakdancer_ssv/$ID.breakdancer.gff"],
			type=>["ssv"],
		};
		$g->add_vertex($_Job);
		$ret{'somatic_sv'}=$_Job;
	}
	return %ret;
}
sub	genReport
{
	my ($pJob)=@_;
	my $report_out=$pJob->{'out'}->[0];
	my $report_in="$opt_o/report/report.in";
	push @{$pJob->{'cmd'}},"$opt_bin/report/report.sh $opt_ini $opt_TR $report_in $report_out";
	mkpath "$opt_o/report/in";

	my %fList=();

	## data production
	`awk -v OFS="\t" '!/^#/{print \$3,\$5,\$9}' $g_infile > $opt_o/report/in/dataProduction.list`;
	$fList{'dataProduction_list'}="$opt_o/report/in/dataProduction.list";
	## stat list
	my $o;
	open $o,">","$opt_o/report/in/stat.list" or die "$!";
	$fList{'stat_list'}="$opt_o/report/in/stat.list";
	foreach (sort keys %{$pJob->{'in'}->{'stat'}})
	{
		print $o "$_\t$pJob->{'in'}->{'stat'}->{$_}->{'out'}->[0]\n";
	}
	## depth list
	open $o,">","$opt_o/report/in/depth.list" or die "$!";
	$fList{'depth_list'}="$opt_o/report/in/depth.list";
	foreach (sort keys %{$pJob->{'in'}->{'stat'}})
	{
		print $o "$_\t$pJob->{'in'}->{'stat'}->{$_}->{'out'}->[1]\n";
	}
	close $o;
	## cumu list
	open $o,">","$opt_o/report/in/cumu.list" or die "$!";
	$fList{'cumu_list'}="$opt_o/report/in/cumu.list";
	foreach (sort keys %{$pJob->{'in'}->{'stat'}})
	{
		print $o "$_\t$pJob->{'in'}->{'stat'}->{$_}->{'out'}->[2]\n";
	}
	close $o;
	## covByChr list
	open $o,">","$opt_o/report/in/covByChr.list" or die "$!";
	$fList{'covByChr_list'}="$opt_o/report/in/covByChr.list";
	foreach (sort keys %{$pJob->{'in'}->{'stat'}})
	{
		print $o "$_\t$pJob->{'in'}->{'stat'}->{$_}->{'out'}->[3]\n";
	}
	close $o;

	## var list
	if($opt_var)
	{
		foreach my $_m (sort keys %{$pJob->{'in'}->{'var'}})
		{
			open $o,">","$opt_o/report/in/var.$_m.list" or die "$!";
			$fList{"var_${_m}_list"}="$opt_o/report/in/var.$_m.list";
			foreach (sort keys %{$pJob->{'in'}->{'var'}->{$_m}})
			{
				print $o "$_\t$pJob->{'in'}->{'var'}->{$_m}->{$_}->{'out'}->[0].stat.txt\n";
				if($pJob->{'in'}->{'var'}->{$_m}->{$_}->{'out'}->[1])
				{
					print $o "$_\t$pJob->{'in'}->{'var'}->{$_m}->{$_}->{'out'}->[1].stat.txt\n";
				}
			}
			close $o;
		}
	}
	if(!defined($opt_single))
	{
		foreach my $_m (sort keys %{$pJob->{'in'}->{'somatic'}})
		{
			open $o,">","$opt_o/report/in/somatic.$_m.list" or die "$!";
			$fList{"somatic_${_m}_list"}="$opt_o/report/in/somatic.$_m.list";
			foreach (sort keys %{$pJob->{'in'}->{'somatic'}->{$_m}})
			{
				print $o "$_\t$pJob->{'in'}->{'somatic'}->{$_m}->{$_}->{'out'}->[0]\n";
				#my $_cmd="echo $pJob->{'in'}->{'somatic'}->{$_m}->{$_}->{'out'}->[0]";
				#push @{$pJob->{'cmd'}},$_cmd;
				if($pJob->{'in'}->{'somatic'}->{$_m}->{$_}->{'out'}->[1])
				{
					print $o "$_\t$pJob->{'in'}->{'somatic'}->{$_m}->{$_}->{'out'}->[1]\n";
					#$_cmd="echo $pJob->{'in'}->{'somatic'}->{$_m}->{$_}->{'out'}->[1]";
					#push @{$pJob->{'cmd'}},$_cmd;
				}
			}
			close $o;
		}
	}
	#printf STDERR "%s\n",join("\n",(@{$pJob->{'cmd'}}));
	open $o,">",$report_in or die "$!";
	foreach (sort keys %fList)
	{
		print $o "$_=$fList{$_}\n";
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
	open $out,">","$_outDir/$_name.sh" or die "$!";
	printf $out "%s\n",join(" && \\\n",@$_cmd);
	close $out;

	## output sjm job
	printf $gOut "job_begin\n";
	printf $gOut "  name $v->{name}\n";
	printf $gOut "  memory $v->{memory}\n";
	if(defined($v->{"cpuNum"})) { printf $gOut "  cpu $v->{cpuNum}\n"; }
	if($opt_queue ne "") { printf $gOut "  queue $opt_queue\n"; }
	if($opt_project ne "") { printf $gOut "  project $opt_project\n"; }
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
