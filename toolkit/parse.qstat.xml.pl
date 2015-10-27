#!/usr/bin/env perl
#============================================================================
# Name        		: parse.qstat.xml.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sun Mar 15 21:54:13 2015
# Last Modified By	: 
# Last Modified On	: Sun Mar 15 21:54:13 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl parse.qstat.xml.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

use lib("/usr/local/lib64/perl5");
use XML::LibXML;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

#print "step1\n";

if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"bgzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}

#print "step2\n";

my $parser=XML::LibXML->new();
my $doc=$parser->parse_file($infile);

print "JB_job_number\tJB_name\tJB_owner\tJB_project\tstat\tNODE_cpu_usage\tNODE_mem_usage\tNODE_io_usage\tNODE_slots\tREQUEST_virtual_free\tREQUEST_num_proc\tREQUEST_queue\tQUEUE_name\tQUEUE_type\tQUEUE_slots_used-resv-total\tNODE_load_avg\tNODE_virtual_free_available\tNODE_num_proc_vailable\n";
foreach my $sge_queue ($doc->findnodes("//Queue-List"))
{
	my $load_avg="";
	my $vf="";
	my $np="";
	foreach my $n ($sge_queue->findnodes("./resource"))
	{
		foreach ($n->attributes())
		{
			if(/name="num_proc"/) 
			{
				$np=$n->textContent();
				last;
			}elsif(/name="virtual_free"/)
			{
				$vf=$n->textContent();
				last;
			}elsif(/name="load_avg"/)
			{
				$load_avg=$n->textContent();
				last;
			}
		}
	}
	foreach my $n ($sge_queue->findnodes("./job_list"))
	{
		my ($request_vf,$request_np)=findRequest($n);
		print 	$n->findnodes("JB_job_number")."\t".
		 	 	$n->findnodes("JB_name")."\t".
		 	 	$n->findnodes("JB_owner")."\t".
		 	 	$n->findnodes("JB_project")."\t".
		 	 	$n->findnodes("state")."\t".
		 	 	$n->findnodes("cpu_usage")."\t".
		 	 	$n->findnodes("mem_usage")."\t".
		 	 	$n->findnodes("io_usage")."\t".
		 	 	$n->findnodes("slots")."\t".
				$request_vf."\t".
				$request_np."\t".
				$n->findnodes("hard_req_queue")."\t".
				$sge_queue->findnodes("./name")."\t".
				$sge_queue->findnodes("./qtype")."\t".
				$sge_queue->findnodes("./slots_used")."/".$sge_queue->findnodes("./slots_resv")."/".$sge_queue->findnodes("./slots_total")."\t".
				$load_avg."\t".
				$vf."\t".
				$np.
				"\n";
	}
}

for my $pending_job ($doc->findnodes("//job_info/job_list"))
{
	my ($request_vf,$request_np)=findRequest($pending_job);
	print $pending_job->findnodes("./JB_job_number")."\t".
			$pending_job->findnodes("./JB_name")."\t".
		 	$pending_job->findnodes("JB_owner")."\t".
		 	$pending_job->findnodes("JB_project")."\t".
		 	$pending_job->findnodes("state")."\t".
			"NA\tNA\tNA\t".
		 	$pending_job->findnodes("slots")."\t".
			$request_vf."\t",
			$request_np."\t",
			$pending_job->findnodes("hard_req_queue")."\t".
			"NA\tNA\tNA\tNA\tNA\tNA".
			"\n";
}

###################################################################

sub findRequest
{
	my ($n)=@_;
	my $request_vf="";
	my $request_np="";
	foreach my $nn ($n->findnodes("./hard_request"))
	{
		foreach ($nn->attributes())
		{
			if(/name="num_proc"/) 
			{
				$request_np=$nn->textContent();
				last;
			}elsif(/name="virtual_free"/)
			{
				$request_vf=$nn->textContent();
				last;
			}
		}
	}
	return ($request_vf,$request_np);
}

sub usage
{
    die `pod2text $0`;
}
