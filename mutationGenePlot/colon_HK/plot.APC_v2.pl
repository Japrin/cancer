#!/usr/bin/perl
#============================================================================
# Name        		: plot.APC.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Thu Nov 11 14:21:18 2010
# Last Modified By	: 
# Last Modified On	: Thu Nov 11 14:21:18 2010
# Copyright   		: Copyright (C) 2010
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl plot.muc17.pl [option] <outfile> <seq file> <ms file> <domain file> <refGene file>

	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use SVG;
use Data::Dumper;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<5 || $opt_h) { usage(); }
my $outfile=shift @ARGV;
my $seqfile=shift @ARGV;
my $msfile=shift @ARGV;
my $domainFile=shift @ARGV;
my $refGeneFile=shift @ARGV;

my $width=650;
my $height=500;
my %aaTable=(
		Ala=>'A', Arg=>'R', Asn=>'N', Asp=>'D', 
		Cys=>'C', Glu=>'E', Gln=>'Q', Gly=>'G',
		His=>'H', Ile=>'I', Leu=>'L', Lys=>'K',
		Met=>'M', Phe=>'F', Pro=>'P', Ser=>'S',
		Thr=>'T', Trp=>'W', Tyr=>'Y', Val=>'V',
		Stop=>'*',
		);

my @lenCDS=();
my @rangeCDS=();
my @posCDS=();
my @domainLen=();
my @domainRegion=();
my @domainPos=();
my @wPerAA=(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,1,0.2,1);
my $autoPerAA=0.1;
my @domainCol=();

my @wPerBase=();
map { push @wPerBase,$_/3; } @wPerAA;
my %glySite=();
printf "======================start:plot.muc17.pl (%s)=================================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);
readRefGene(\@lenCDS,\@rangeCDS,\@posCDS,$refGeneFile);
readDomain(\@domainLen,\@domainRegion,\@domainPos,$domainFile);
#print Dumper @domainCol;
my $canvas = SVG->new('width',$width,'height',$height);
$canvas->rect(x=>0,y=>0,width=>$width,height=>$height,fill=>"white");

my $m=$canvas->marker('id','arrow','viewBox',"0 0 10 10",'refX',"0",'refY',"5",'markerUnits',"strokeWidth",'markerWidth',"8",'markerHeight',"6",'orient',"auto");
$m->path('d',"M 0 0 L 15 5 L 0 10 z");
my $m2=$canvas->marker('id','arrow2','viewBox',"0 0 10 10",'refX',"0",'refY',"5",'markerUnits',"strokeWidth",'markerWidth',"8",'markerHeight',"6",'orient',"auto",'transform',"rotate(180)");
$m2->path('d',"M 15 0 L 0 5 L 15 10 z");

my $gfill1=$canvas->tag('radialGradient','id','myRadialGradient1','r','65%');
$gfill1->tag('stop','offset','5%','stop-color','red');
$gfill1->tag('stop','offset','95%','stop-color','darkred');
my $gfill2=$canvas->tag('linearGradient','id','myLinearGradient2','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
$gfill2->tag('stop','offset','5%','stop-color','blue');
$gfill2->tag('stop','offset','50%','stop-color','lightskyblue');
$gfill2->tag('stop','offset','95%','stop-color','blue');
my $gfill3=$canvas->tag('linearGradient','id','myLinearGradient3','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
$gfill3->tag('stop','offset','5%','stop-color','green');
$gfill3->tag('stop','offset','50%','stop-color','lightgreen');
$gfill3->tag('stop','offset','95%','stop-color','green');
my $gfill4=$canvas->tag('linearGradient','id','myLinearGradient4','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
$gfill4->tag('stop','offset','5%','stop-color','deeppink');
$gfill4->tag('stop','offset','50%','stop-color','lightpink');
$gfill4->tag('stop','offset','95%','stop-color','deeppink');

###
#star command
my $starCol="green";
my $starCmd="M -0.951,-0.309 L  0.951,-0.309 -0.588, 0.809 0.000,-1.000 0.588, 0.809 Z";

my $margin=20;
plotGeneStruct($canvas,$margin,250);
plotSNVBySeq($seqfile,$canvas,$margin,250);
plotSNVByMS($msfile,$canvas,$margin,280);
plotLegend($canvas,440,30);
open $out,">",$outfile or die "Cann't open file $outfile ($!)\n";
print $out $canvas->xmlify;

printf "======================finished:plot.muc17.pl (%s)==============================\n",strftime("%a %b %e %H:%M:%S %Y", localtime);



############################################################################
sub readDomain
{
	#readDomain(\@domainLen,\@domainRegion,\@domainPos,$domainFile);
	my ($listLen,$listRange,$listPos,$infile)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	my $cumLen=0;
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @field=split /\t/;
		#Q685J3  UniProtKB       Signal peptide  1       25      .       .       .       Status=Potential
		if(@$listRange && $listRange->[-1]+1<$field[3])
		{
			push @$listLen,$field[3]-$listRange->[-1]-1;
			push @$listRange,$listRange->[-1]+1,$field[3]-1;
			push @$listPos,$cumLen+1;
			$cumLen+=$field[3]-$listRange->[-1]-1;
			push @domainCol,"rgb(97,199,199)";
		}
		push @$listLen,$field[4]-$field[3]+1;
		push @$listRange,$field[3],$field[4];
		push @$listPos,$cumLen+1;
		$cumLen+=$field[4]-$field[3]+1;
		if($field[2] eq 'Coiled coil') { push @domainCol,"hotpink"; }
		elsif($field[2] eq 'Repeat') { push @domainCol,"lightgreen"; }
		elsif($field[2] eq 'Region') { push @domainCol,'rgb(233,150,122)'; }
		elsif($field[2] eq 'Motif') { push @domainCol,'red'; }
		else { warn "$_\n"; }
	}
}
sub plotSNVByMS
{
	my ($infile,$canvas,$x_offset,$y)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	my @dd=();
	my @ii=();
	my @nn=();
	my @ff=();
	my @gg=();
	while(<$in>)
	{
		chomp;
		if(/^\s*$/) { next; }
		my @field=split ;
		#chr7	100461207	TAG=>AAG	57;	Stop=>Lys	readthrough	110NT;
		my $basePos=$field[1];
		$field[3]=~/(.+);/;
		my $aaPos=$1;
		if($aaPos eq 'NA' || $aaPos == -1 ) { $aaPos=aaPos($basePos,\@posCDS,\@rangeCDS); }
		my @tt=split(/;/,$field[6]);
		my $nSample=@tt;
		my $mutType=$field[5];
		##
		my $pos=$x_offset;
		my $i=whichRegion($aaPos,\@domainRegion);
		for(my $j=0;$j<$i-1;$j++)
		{
			$pos+=$wPerAA[$j]*$domainLen[$j];
		}
		$pos+=$wPerAA[$i-1]*($aaPos-$domainRegion[($i-1)*2]);
		##
		$canvas->line('x1',$pos,'y1',$y,'x2',$pos,'y2',$y+20,'style',{'stroke','black'});
		my $fillCol=($mutType eq 'nonsense'?'red':'mediumslateblue');
		my $text=$field[4];
		$text=~/(.+)=>(.+)/;
		$text=sprintf("%s%s%s",$aaTable{$1},$aaPos,$aaTable{$2});

		push @ii,$pos;
		push @nn,$nSample;
		push @ff,$fillCol;
		push @dd,$text;
		push @gg,$basePos;
	}
	my $txt_x=41;
	my $txt_y=487;
	for(my $i=0;$i<@dd;$i++)
	{
		my $text=$dd[$i];
		my $pp=$txt_x-5;
		if($i==0) { $pp=$ii[$i]-15; }
		elsif($i==1) { $pp=$ii[$i]; }
		elsif($i==2) { $pp=$ii[$i]+10; }
		elsif($i==3) { $pp=$ii[$i]; }
		elsif($i==4) { $pp=$ii[$i]-13; }
		elsif($i==5) { $pp=$ii[$i]; }
		elsif($i==6) { $pp=$ii[$i]+13; }
		elsif($i==7) { $pp=$ii[$i]; }
		$txt_x=$pp+2;
		my $txtColor="black";
		$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','12','-cdata',"$text",'fill',"$txtColor",'transform',"translate($txt_x,$txt_y) rotate(-90)");
		#$txt_x+=19;
		my $cmd=sprintf("M %s %s L %s %s L %s %s ",$ii[$i],$y+20,$pp,$y+40,$pp,$y+150);
		$canvas->path('d',$cmd,'style',{'stroke','black','fill','none','marker-end',"url(#arrow)"});
		for(my $k=0;$k<$nn[$i];$k++)
		{
			$canvas->circle('cx',$pp,'cy',$y+50+$k*7,'r',3,'style',{'fill',$ff[$i]});
		}
		#print $gg[$i],"\n";
		if(exists($glySite{$gg[$i]}))
		{
			$canvas->path('d',$starCmd,'style',{'fill',$starCol},'transform',"translate($pp,448) scale(6,6)");
		}
	}
}
sub plotSNVBySeq
{
	my ($infile,$canvas,$x_offset,$y)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	my @dd=();
	my @ii=();
	my @nn=();
	my @ff=();
	my @gg=();
	while(<$in>)
	{
		chomp;
		if(/^\s*$/) { next; }
		my @field=split ;
		#chr7	100461207	TAG=>AAG	57;	Stop=>Lys	readthrough	110NT;
		my $basePos=$field[1];
		$field[3]=~/(.+);/;
		my $aaPos=$1;
		if( $aaPos eq 'NA' || $aaPos == -1 ) { $aaPos=aaPos($basePos,\@posCDS,\@rangeCDS); }
		my @tt=split(/;/,$field[6]);
		my $nSample=@tt;
		my $mutType=$field[5];
		##
		my $pos=$x_offset;
		my $i=whichRegion($aaPos,\@domainRegion);
		for(my $j=0;$j<$i-1;$j++)
		{
			$pos+=$wPerAA[$j]*$domainLen[$j];
		}
		$pos+=$wPerAA[$i-1]*($aaPos-$domainRegion[($i-1)*2]);
		##
		$canvas->line('x1',$pos,'y1',$y,'x2',$pos,'y2',$y-20,'style',{'stroke','black'});
		my $fillCol=($mutType eq 'nonsense'?'red':'mediumslateblue');
		my $text=$field[4];
		$text=~/(.+)=>(.+)/;
		$text=sprintf("%s%s%s",$aaTable{$1},$aaPos,$aaTable{$2});

		push @ii,$pos;
		push @nn,$nSample;
		push @ff,$fillCol;
		push @dd,$text;
		push @gg,$basePos;
	}
	my $txt_x=41;
	my $txt_y=88;
	for(my $i=0;$i<@dd;$i++)
	{
		my $text=$dd[$i];
		my $pp=$txt_x-5;
		if($i==0) { $pp=$ii[$i]; }
		elsif($i==1) { $pp=$ii[$i]; }
		elsif($i==2) { $pp=$ii[$i]+10; }
		elsif($i==3) { $pp=$ii[$i]; }
		elsif($i==4) { $pp=$ii[$i]-10; }
		elsif($i==5) { $pp=$ii[$i]+10; }
		elsif($i==6) { $pp=$ii[$i]; }
		elsif($i==7) { $pp=$ii[$i]; }
		$txt_x=$pp+2;
		my $txtColor="black";
		$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','12','-cdata',"$text",'fill',"$txtColor",'transform',"translate($txt_x,$txt_y) rotate(-90)");
		#$txt_x+=19;
		my $cmd=sprintf("M %s %s L %s %s L %s %s ",$ii[$i],$y-20,$pp,$y-40,$pp,$y-150);
		$canvas->path('d',$cmd,'style',{'stroke','black','fill','none','marker-end',"url(#arrow)"});
		for(my $k=0;$k<$nn[$i];$k++)
		{
			$canvas->circle('cx',$pp,'cy',$y-50-$k*7,'r',3,'style',{'fill',$ff[$i]});
		}
		#print $gg[$i],"\n";
		if(exists($glySite{$gg[$i]}))
		{
			$canvas->path('d',$starCmd,'style',{'fill',$starCol},'transform',"translate($pp,448) scale(6,6)");
		}
	}
}
sub plotLegend
{
	my ($canvas,$x,$y)=@_;
	$canvas->rect('x',$x,'y',$y,'width',150,'height',160,'style',{'fill','white','stroke','black','stroke-width',1.5});
	$canvas->circle('cx',$x+30,'cy',$y+11,'r',3,'style',{'fill','red'});
	$canvas->text('x',$x+60,'y',$y+14,'font-family','TimesNewRoman','font-size','12','-cdata',"Nonsense",'fill','black');
	$canvas->circle('cx',$x+30,'cy',$y+26,'r',3,'style',{'fill','mediumslateblue'});
	$canvas->text('x',$x+60,'y',$y+29,'font-family','TimesNewRoman','font-size','12','-cdata',"Missense",'fill','black');

	$canvas->rect('x',$x+10,'y',$y+35,'width',40,'height',25,'style',{'fill','hotpink','stroke','black'});
	$canvas->text('x',$x+60,'y',$y+50,'font-family','TimesNewRoman','font-size','12','-cdata',"Coiled coil",'fill','black');
	$canvas->rect('x',$x+10,'y',$y+65,'width',40,'height',25,'style',{'fill','lightgreen','stroke','black'});
	$canvas->text('x',$x+60,'y',$y+80,'font-family','TimesNewRoman','font-size','12','-cdata',"Repeat",'fill','black');
	$canvas->rect('x',$x+10,'y',$y+95,'width',40,'height',25,'style',{'fill','rgb(233,150,122)','stroke','black'});
	$canvas->text('x',$x+60,'y',$y+110,'font-family','TimesNewRoman','font-size','12','-cdata',"Region[1]",'fill','black');
	$canvas->rect('x',$x+10,'y',$y+125,'width',40,'height',25,'style',{'fill','red','stroke','black'});
	$canvas->text('x',$x+60,'y',$y+140,'font-family','TimesNewRoman','font-size','12','-cdata',"Motif[2]",'fill','black');
	
	#my $star_x=$x+15;
	#my $star_y=$y+41;
	#$canvas->path('d',$starCmd,'style',{'fill',$starCol},'transform',"translate($star_x,$star_y) scale(6,6)");
	#$canvas->text('x',$x+30,'y',$y+44,'font-family','TimesNewRoman','font-size','12','-cdata',"glycosylation site",'fill','black');
}
sub plotGeneStruct
{
	my ($canvas,$x,$y)=@_;
	$height=30;
	my $p=$x;
	my $w=0;
	for (my $i=0;$i<@domainLen;$i++)
	{
		$w=$wPerAA[$i]*$domainLen[$i];
		print "$wPerAA[$i]\t$domainLen[$i]\t$w\n";
		my $_domainCol=$domainCol[$i];
		$canvas->rect('x',$p,'y',$y,'width',$w,'height',$height,'style',{'fill',"$_domainCol",'stroke','black'});
		$p+=$w;
	}
	my @scale=(0,400,800,1200,1600,2000,2400,2800);
	foreach my $_x (@scale)
	{
		my $d_x=$x;
		for(my $_i=0;$_i<@domainPos;$_i++)
		{
			if($domainPos[$_i]<=$_x || $_x==0)
			{
				$d_x+=$wPerAA[$_i]*($_x-$domainPos[$_i]+1);
				$canvas->line('x1',$d_x,'y1',$y+$height,'x2',$d_x,'y2',$y+$height+12,'style',{'stroke','black'});
				$canvas->text('x',$d_x,'y',$y+$height+24,'font-family','TimesNewRoman','font-size','14','-cdata',$_x,'fill','black','stroke-width',1.5,'text-anchor','middle');
				last;
			}else
			{
				$d_x+=$wPerAA[$_i]*($domainLen[$_i]);
			}
		}
	}
	#
}
sub whichRegion
{
	my ($aaPos,$listRange)=@_;
	for(my $i=0;$i<(@$listRange)/2;$i++)
	{
		my $pos1=$listRange->[$i*2];
		my $pos2=$listRange->[$i*2+1];
		if($pos1<=$aaPos && $aaPos<=$pos2)
		{
			return $i+1;
		}
	}
	my $str=sprintf("$aaPos no in region:[%s]\n",join(",",@$listRange));
	warn $str;
	return 0;
}
sub aaPos
{
	my ($basePos,$listPos,$listRange)=@_;
	for(my $i=0;$i<@$listPos;$i++)
	{
		my $pos1=$listRange->[$i*2];
		my $pos2=$listRange->[$i*2+1];
		if($pos1<=$basePos && $basePos<=$pos2)
		{
			my $ret=int(($basePos-$pos1+$listPos->[$i]-1)/3)+1;
			#print "basePos:$basePos\taaPos:$ret\n";
			return $ret;
		}
	}
	#my $str=sprintf("$basePos no in CDS:[%s]\n",join(",",@$listRange));
	#warn $str;
	return 0;
}
sub readRefGene
{
	my ($listLen,$listRange,$listPos,$infile)=@_;
	my $in;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	my $cumLen=0;
	while(<$in>)
	{
		chomp;
		if(/^\s*$/) { next; }
		my @field=split;
		#chr7    refSeq  CDS     100450137       100450218       .       +       0       Parent=NM_001040105;
		if($field[2]!~/CDS/) { next; }
		#$list
		push @$listLen,$field[4]-$field[3]+1;
		push @$listRange,$field[3],$field[4];
		push @$listPos,$cumLen+1;
		$cumLen+=$field[4]-$field[3]+1;
	}
}
sub usage
{
	die `pod2text $0`;
}
