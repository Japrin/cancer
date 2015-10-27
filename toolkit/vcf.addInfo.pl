#!/usr/bin/perl
#============================================================================
# Name        		: vcf.addInfo.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Thu Sep 12 09:58:11 2013
# Last Modified By	: 
# Last Modified On	: Thu Sep 12 09:58:11 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl vcf.addInfo.pl [option] [infile]
	
	-l	list file
	-j	key column in infile, "," seperated string [default "0"]
	-i	key column in list file, "," seperated string [default "0"]
	-v	value column in list file, 0~$#F. [default -1,no]
	-u	external value [used only when -v is -1]
	-s	value lable [default "myInfo"]
	-t	value column in infile [default 7]
	-d	maximun distance between sites allowed for sites to be considered overlap [default 0]
	-a	auto extend, so -d will be non't effective [default OFF]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_l,$opt_i,$opt_v,$opt_u,$opt_s,$opt_d,$opt_j,$opt_t,$opt_a);
GetOptions("h"	=>\$opt_h,"l=s"=>\$opt_l,"i=s"=>\$opt_i,"v=i"=>\$opt_v,"u=s"=>\$opt_u,"s=s"=>\$opt_s,"d=i"=>\$opt_d,"j=s"=>\$opt_j,"t=i"=>\$opt_t,"a"=>\$opt_a);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_l)) { usage(); }
if(!defined($opt_i)) { $opt_i=0; }
if(!defined($opt_j)) { $opt_j="0,1"; }
if(!defined($opt_v)) { $opt_v=-1; }
if(!defined($opt_u)) { $opt_u=""; }
if(!defined($opt_s)) { $opt_s="myInfo"; }
if(!defined($opt_d)) { $opt_d=0; }
if(!defined($opt_t)) { $opt_t=7; }

my @index_i=split /,/,$opt_i;
my @index_j=split /,/,$opt_j;
#print STDERR "opt_i\t$opt_i\n";
#exit 0;

my %list=();
readList(\%list,$opt_l);
while(<>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { print "$line\n"; next; }
	my @F=split /\t/;
	my ($chr,$pos)=@F[@index_j];
	
	my $ref=$F[$index_j[1]+2];
	my $alt=$F[$index_j[1]+3];
	my $len_m=length($ref);
	my $len_n=length($alt);
	my $ll=$len_m>$len_n?$len_m:$len_n;
	my ($_l,$_r);
	if($opt_a)
	{
		$_l=$pos;
		$_r=$pos+$ll-1;
	}else
	{
		$_l=$pos-$opt_d;
		$_r=$pos+$opt_d;
	}
	
	for(my $_p=$_l;$_p<=$_r;$_p++)
	{
		if($list{"$chr:$_p"})
		{
			if($F[$opt_t] eq ".") { $F[$opt_t]="$opt_s=".$list{"$chr:$_p"}; }
			else {  $F[$opt_t].=";$opt_s=".$list{"$chr:$_p"};}
			last;
		}
	}
	print join("\t",@F)."\n";
}
###################################################################



sub readList
{
	my $in;
	my ($pList,$infile)=@_;
	if($infile=~/\.gz$/) { open $in,"gzip -cd $infile | " or die "$!"; }
	else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
	while(<$in>)
	{
		chomp;
		my $line=$_;
		if(/^\s*$/ || /^#/) { next; }
		my @F=split /\t/;
		#print STDERR join("\t",@index_i)."\n";
		$F[$index_i[0]]=~s/^chr//;
		my $key=join(":",@F[@index_i]);
		my $v;
		if($opt_v == -1)
		{
			$v=$opt_u;
		}else
		{
			$v=$F[$opt_v];
		}
		$pList->{$key}=$v;
	}
}

sub usage
{
	die `pod2text $0`;
}
