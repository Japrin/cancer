#!/usr/bin/env perl
#============================================================================
# Name        		: venn7.pl
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

	perl venn7.pl [option] <outDir> <infile1> <infile2> <infile3>
	
	-i	id in <infiles> [default 0] (can be a string seprated by comma)
	-n1=<str>
	-n2=<str>
	-n3=<str>
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
my ($opt_h,$opt_i,$opt_j,$opt_n1,$opt_n2,$opt_n3);
GetOptions("h"	=>\$opt_h,"i=s"=>\$opt_i,"n1=s"=>\$opt_n1,"n2=s"=>\$opt_n2,"n3=s"=>\$opt_n3);

if(@ARGV<4 || $opt_h) { usage(); }
my $outDir=shift @ARGV;
my $infile1=shift @ARGV;
my $infile2=shift @ARGV;
my $infile3=shift @ARGV;
$outDir=abs_path($outDir);

if(!defined($opt_i)) { $opt_i="0"; }
$opt_n1 ||=basename($infile1);
$opt_n2 ||=basename($infile2);
$opt_n3 ||=basename($infile3);

printf "======================start:venn7.pl (%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
my %list=();
my @arg_i=split(/,/,$opt_i);
readList(\%list,$infile1,1,@arg_i);
readList(\%list,$infile2,2,@arg_i);
readList(\%list,$infile3,4,@arg_i);
#print Dumper(%list);exit(0);
my %part=(1=>[],2=>[],3=>[],4=>[],5=>[],6=>[],7=>[]);
foreach my $key (sort keys %list)
{
	my $ii=$list{$key}->{'part'};
	my $ss=join(" | ",@{$list{$key}->{'entry'}});
	push @{$part{$ii}},"$key\t$ss";
}
my @lable=();
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
plotSVG(\@lable,$opt_n1,$opt_n2,$opt_n3,$outDir);
printf "======================finished:venn7.pl (%s)==============================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);



############################################################################
sub plotSVG
{
	my ($pLable,$n1,$n2,$n3,$outDir)=@_;
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
	#draw the circle of $file3
	$venn -> circle( 
		cx => $cx_1+$dist/2,
		cy => $cy_1+sqrt(3*($dist*$dist)/4),
		r => $r,
		fill => 'green',
		stroke => 'green',
		'stroke-width' => $stroke_width,
		opacity => $opacity,
	);
	my $text_fill = 'black';
	my $font_size = int($r/5);
	my $text_anchor = 'middle';
	my $move = 20;
	$venn -> text( # write the number of part 1
		x => $cx_1-$move,
		y => $cy_1,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[0]);
	$venn -> text( # write the number of part 2
		x => $cx_1+$dist+$move,
		y => $cy_1,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[1]);
	$venn -> text( # write the number of part 3
		x => $cx_1+$dist/2,
		y => $cy_1-$move,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[2]);
	$venn -> text( # write the number of part 4
		x => $cx_1+$dist/2,
		y => $cy_1+sqrt(3*($dist*$dist)/4)+$move,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[3]);
	$venn -> text( # write the number of part 5
		x => $cx_1,
		y => $cy_1+2*$r/3,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[4]);
	$venn -> text( # write the number of part 6
		x => $cx_1+$dist,
		y => $cy_1+2*$r/3,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[5]);
	$venn -> text( # write the number of part 7
		x => $cx_1+$dist/2,
		y => $cy_1+$r/2,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($pLable->[6]);
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
	$venn -> text( # write the name of circle3
		x => $cx_1+$dist/2,
		y => $cy_1+sqrt(3*($dist*$dist)/4)+$r+$move,
		fill => $text_fill,
		'font-size' => $font_size,
		'text-anchor' => $text_anchor,
	) -> cdata($n3);
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
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @field=split /\t/;
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
