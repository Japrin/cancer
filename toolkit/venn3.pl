#!/usr/bin/env perl
#============================================================================
# Name        		: venn3.pl
# Author      		: zhenglt
# Version     		: v1.00
# Created On  		: Sun Jan 16 14:15:53 2011
# Last Modified By	: 
# Last Modified On	: Sun Jan 16 14:15:53 2011
# Copyright   		: Copyright (C) 2011
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl venn3.pl [option] <outDir> [<infile1> <infile2> | <num1> <num2> <num3> ]
	
	-i	id in <infiles> [default 0] (can be a string seprated by comma)
	-n1	<str> name of circle1
	-n2	<str> name of circle2
	-title	<str> title
	-sep	seperator [default {tab}]
	-header	list with tab-header
	-t	reads number from num1, num2, num3 (not from input file)
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;
use File::Basename;
use Cwd qw(abs_path);
use Data::Dumper;
use SVG;

my ($in,$out);
my ($opt_h,$opt_i,$opt_j,$opt_n1,$opt_n2,$opt_sep,$opt_header,$opt_t,$opt_title);
GetOptions("h"	=>\$opt_h,
			"i=s"=>\$opt_i,
			"n1=s"=>\$opt_n1,
			"n2=s"=>\$opt_n2,
			"sep=s"=>\$opt_sep,
			"header"=>\$opt_header,
			"t"		=>\$opt_t,
			"title=s"	=>\$opt_title,
			);

if(@ARGV<3 || $opt_h) { usage(); }
my $outDir=shift @ARGV;
mkpath $outDir;
$outDir=abs_path($outDir);

if(!defined($opt_title)) { $opt_title="Venn figure"; }

printf "======================start:venn3.pl (%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
my %list=();
my @lable=();
if($opt_t)
{
	@lable=@ARGV;
}else
{
	my $infile1=shift @ARGV;
	my $infile2=shift @ARGV;
	if(!defined($opt_n1)) { $opt_n1=basename($infile1); }
	if(!defined($opt_n2)) { $opt_n2=basename($infile2); }
	if(!defined($opt_i)) { $opt_i="0"; }
	if(!defined($opt_sep)) { $opt_sep="\t"; }

	my @arg_i=split(/,/,$opt_i);
	readList(\%list,$infile1,1,@arg_i);
	readList(\%list,$infile2,2,@arg_i);
#print Dumper(%list);exit(0);
	my %part=(1=>[],2=>[],3=>[]);
	foreach my $key (sort keys %list)
	{
		my $ii=$list{$key}->{'part'};
		my $ss=join(" | ",@{$list{$key}->{'entry'}});
		push @{$part{$ii}},"$key\t$ss";
	}
	mkpath $outDir;
	foreach (sort { $a<=>$b } keys %part)
	{
		open $out,">","$outDir/part$_" or die "Cann't open file $outDir/part$_ ($!) \n";
		printf ">part$_: %d\n",scalar @{$part{$_}};
		push @lable,scalar @{$part{$_}};
		foreach my $s (@{$part{$_}})
		{
			print $out "$s\n";
		}
	}
}
plotSVG(\@lable,$opt_n1,$opt_n2,$opt_title,$outDir);
printf "======================finished:venn3.pl (%s)==============================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);



############################################################################
sub plotSVG
{
	my ($pLable,$n1,$n2,$title,$outDir)=@_;
	my ($width, $height) = (640, 400);
	my $venn = SVG -> new( width => $width, height => $height,);
	my $r = 100;
	my $dist = 1.2*$r;
	my $cx_1 = ($width-$dist)/2;
	my $cy_1 = ($height-sqrt(3*($dist*$dist)/4))/2;
	my ($opacity, $stroke_width) = (0.5, 1);
	$venn->rect('x',0,'y',0,'width',$width,'height',$height,'fill','white');
	#draw the circle of $file1
	$venn -> circle(
		cx => $cx_1,
		cy => $cy_1,
		r => $r,
		fill => 'red',
		stroke => 'red',
		'stroke-width' => $stroke_width,
		opacity => $opacity,
	);
	#draw the circle of $file2
	$venn -> circle(
		cx => $cx_1+$dist,
		cy => $cy_1,
		r => $r,
		fill => 'blue',
		stroke => 'blue',
		'stroke-width' => $stroke_width,
		opacity => $opacity,
	);
	my $text_fill = 'black';
	my $font_size = int($r/5)-2;
	my $text_anchor = 'middle';
	my $move = 20;
	my @_t=split /[\(\)]/,$pLable->[0];
	$venn -> text( # write the number of part 1
		x => $cx_1-$move-20,
		y => $cy_1,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($_t[0]);
	if(defined($_t[1])) { $venn->text(x=>$cx_1-$move/2-30,y=>$cy_1+20,fill=>$text_fill,'font-size' => $font_size,'text-anchor' => $text_anchor,)->cdata("($_t[1])"); }
	@_t=split /[\(\)]/,$pLable->[1];
	$venn -> text( # write the number of part 2
		x => $cx_1+$dist+$move+20,
		y => $cy_1,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($_t[0]);
	if(defined($_t[1])) { $venn->text(x=>$cx_1+$dist+$move+20,y=>$cy_1+20,fill=>$text_fill,'font-size' => $font_size,'text-anchor' => $text_anchor,)->cdata("($_t[1])"); }
	@_t=split /[\(\)]/,$pLable->[2];
	$venn -> text( # write the number of part 3
		x => $cx_1+$dist/2,
		y => $cy_1,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($_t[0]);
	if(defined($_t[1])) { $venn->text(x=>$cx_1+$dist/2,y=>$cy_1+20,fill=>$text_fill,'font-size' => $font_size,'text-anchor' => $text_anchor,)->cdata("($_t[1])"); }
##############################################
	$venn -> text( # write the name of circle1
		x => $cx_1-3*$move,
		y => $cy_1-$r-$move/2,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($n1);
	$venn -> text( # write the name of circle2
		x => $cx_1+$dist+3*$move,
		y => $cy_1-$r-$move/2,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($n2);
	$venn -> text( # write the title
		x => $width/2,
		y => 300,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($title);
	open PLOT,">$outDir/venn.svg";
	print PLOT $venn -> xmlify,"\n";
	close PLOT;
}

sub readList
{
	my $list=shift @_;
	my $infile=shift @_;
	my $n=shift @_;
	my @index=@_;
	my $in;
	if($infile=~/\.vcf\.gz$/) { open $in,"bgzip -cd $infile | " or die "$!"; }
	elsif($infile=~/\.gz$/) { open $in,"gzip -cd $infile | " or die "$!"; }
	else { open $in,$infile or die "Cann't open file $infile ($!) \n"; }
	if($opt_header) { my $_temp=<$in>; }
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @field=split /$opt_sep/;
		my $malformed=0;
		for(my $i=0;$i<@index;$i++) { if(!defined($field[$index[$i]])) { $malformed=1; } }
		if($malformed) { next; }
		my $key=join(":",@field[@index]);
		#print "$key\n";
		#my $key=(split)[$index];
		if(!exists($list->{$key}))
		{
			$list->{$key}={part=>0,entry=>[]};
		}
		$list->{$key}->{'part'}+=$n;
		push @{$list->{$key}->{'entry'}},$_;
	}
}
sub usage
{
	die `pod2text $0`;
}
