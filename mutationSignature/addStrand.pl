#!/usr/bin/perl
#============================================================================
# Name        		: addStrand.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Mon Apr 22 18:02:04 2013
# Last Modified By	: 
# Last Modified On	: Mon Apr 22 18:02:04 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl addStrand.pl [option] [infile]

	-l	strand file [required]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_l);
GetOptions("h"	=>\$opt_h,"l=s"=>\$opt_l);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_l)) { usage(); }

my %list=();
readList(\%list,$opt_l);
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { next; }
	my @F=split /\t/;
	#chr2    25834820        .       C       G       .       PASS    Func=intronic;Gene=DTNB;DP=67;SOMATIC;SS=2;SSC=20;GPV=1E0;SPV=8.6346E-3
	#chr2    26712125        .       C       G       .       PASS    Func=exonic;Gene=OTOF;ExonicFunc=synonymousSNV;AAChange=OTOF:NM_194248:exon11:c.G999C:p.L333L
	if($F[7]=~/Func=intergenic;/) { next; }	
	my ($chr,$pos,$ref,$alt)=@F[0,1,3,4];
	my ($gname) = $F[7]=~/Gene=(.+?);/;
	$gname=(split /[,\(]/,$gname)[0];
	if(!exists($list{"${chr}:${gname}"})) { printf STDERR "no exists in strand file: $chr, $gname\n";next; }
	if($list{"${chr}:${gname}"}->{'uniq'}==1)
	{
		printf "$chr\t$pos\t$ref\t$alt\t$gname\t%s\n",$list{"${chr}:${gname}"}->{'strand'}->[0];
	}else
	{
		my $_p=$list{"${chr}:${gname}"}->{'beg'};
		my $_q=$list{"${chr}:${gname}"}->{'end'};
		for(my $i=0;$i<@$_p;$i++)
		{
			if($_p->[$i]<=$pos && $pos<=$_q->[$i])
			{
				printf "$chr\t$pos\t$ref\t$alt\t$gname\t%s\n",$list{"${chr}:${gname}"}->{'strand'}->[$i];
				last;
			}
		}
	}
}
###################################################################



sub readList
{
	my $in;
	my ($pList,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		my ($chr,$beg,$end,$gname,$tname,$strand)=@F;
		#chr1    1115076 1121243 TTLL10  NM_153254       +
		#chr1    1138887 1142089 TNFRSF18        NM_004195       -
		if(!exists($pList->{"${chr}:${gname}"})) { $pList->{"${chr}:${gname}"}={'strand'=>[],'transcript'=>[],'c'=>[0,0],'beg'=>[],'end'=>[],'uniq'=>-1}; };
		push @{$pList->{"${chr}:${gname}"}->{'strand'}},$strand;
		push @{$pList->{"${chr}:${gname}"}->{'transcript'}},$tname;
		push @{$pList->{"${chr}:${gname}"}->{'beg'}},$beg;
		push @{$pList->{"${chr}:${gname}"}->{'end'}},$end;
		if($strand eq "+") { $pList->{"${chr}:${gname}"}->{'c'}->[0]++; } 
		else { $pList->{"${chr}:${gname}"}->{'c'}->[1]++; }
	}
	foreach my $g (keys %$pList)
	{
		if($pList->{$g}->{'c'}->[0] != @{$pList->{$g}->{'strand'}} && $pList->{$g}->{'c'}->[1] != @{$pList->{$g}->{'strand'}})
		{
			$pList->{$g}->{'uniq'}=0;
			#printf STDERR "$g\n";
		}else
		{
			$pList->{$g}->{'uniq'}=1;
		}
	}
}

sub usage
{
	die `pod2text $0`;
}
