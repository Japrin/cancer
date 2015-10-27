#!/PROJ/GR/share/Software/perl/bin/perl
#============================================================================
# Name        		: svg_SV.pl
# Author      		: kongguanyi & zhengliangtao
# Version     		: v1.00
# Created On  		: Wed Nov 13 14:29:58 2013
# Last Modified By	: 
# Last Modified On	: Wed Nov 13 14:29:58 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl svg_SV.pl [option] 
	-i	<region bed> 
	-o	<output prefix>
	-b	<breakdancer>
	-c	<cnv>
	-m	<mappability>
	-p	<2 breakpoints to hightlight>
	-gc	<ref> [default "/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/human_g1k_v37_decoy.fasta"]
	-g	<gene> [default "/PROJ/GR/share/medinfo.00database/annovar/humandb_b37/hg19_refGene.txt"]
	-e	<mappability encoding> [default "/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/gem/human_g1k_v37_decoy.mappability.100bp.out.encoding"]
	-d	<average depth> [default 30]
	-h	display this help and exit

=cut


use strict;
use warnings;
use Getopt::Long;
use SVG;
use Bio::DB::Sam;

my ($opt_i,$out,$out_prefix,$break,$cnv,$posi,$genome,$gene,$opt_h,$opt_m,$opt_e,$opt_RD);
GetOptions(
  "i|in-file=s"  =>\$opt_i,
  "o|out-file=s"    =>\$out_prefix,
  "b=s" =>\$break,
  "c=s" =>\$cnv,
  "p=s" =>\$posi,
  "gc=s" =>\$genome,
  "g=s" =>\$gene,
  "h"   =>\$opt_h,
  "m=s"	=>\$opt_m,
  "e=s"	=>\$opt_e,
  "d=i" =>\$opt_RD,
  );

if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($genome)) { $genome="/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/human_g1k_v37_decoy.fasta"; }
if(!defined($gene)) { $gene="/PROJ/GR/share/medinfo.00database/annovar/humandb_b37/hg19_refGene.txt"; }
if(!defined($opt_e)) { $opt_e="/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/gem/human_g1k_v37_decoy.mappability.100bp.out.encoding"; }
if(!defined($opt_RD)) { $opt_RD=30; }
 
my $fai = Bio::DB::Sam::Fai->load($genome);

my %ENCODING=();
readEncoding(\%ENCODING,$opt_e);

#region to plot
open IN, $opt_i or die "ERROR: open $opt_i: $!";
my @posi=split(/,/,$posi);

my $g_i=0;

while(defined(my $line=<IN>))
{
	chomp $line;
	$g_i++;
	my @colu=split(/\s+/,$line);
	my $chr=$colu[0];
	$chr=~s/^chr//;
	my $start=$colu[1];
	my $end=$colu[2];
	
	print ">begin: $chr:$start-$end\n";

	## g_step: window size in genome; g_resolution: how many pixes in svg per window
	my $g_step=100;
	my $g_resolution=1;

	my $g_NWin=int(($end-$start)/$g_step);
	my $g_width=800;
	my $g_height=800;
	if($g_NWin>2000)
	{
		$g_step=int(($end-$start)/2000);
		my $_aa=$g_step%100;
		$g_step=$g_step-$_aa;
		$g_NWin=int(($end-$start)/$g_step);
		$g_width=$g_NWin+200;
	}elsif($g_NWin<200)
	{
		$g_resolution=10;
	}
	print ">number of plot window: $g_NWin\n";
	print ">window size: $g_step\n";

	my ($BREAK,$CNV,$GENE,$MAP);
	open $BREAK, $break or die "ERROR: open $!";
	open $CNV, $cnv or die "ERROR: open $!";
	open $GENE, $gene or die "ERROR: open $!";
	open $MAP, $opt_m or die "ERROR: open $!";

	### SVG 
	print ">width: $g_width\theight: $g_height\n";
	open $out, ">$out_prefix.$chr.$start.$end.svg" or die "ERROR: open >$out_prefix.$chr.$start.$end: $!";
	my $svg = SVG->new('width',$g_width,'height',$g_height);
	
	## PEM plot
	print ">PEM\n";
	plotPEM($svg,100,150,$BREAK,$chr,$start,$end,$g_step,$g_resolution);

	### RD
	print ">RD\n";
	plotRD($svg,100,300,$CNV,$chr,$start,$end,$g_step,$g_resolution);

	## GC content
	print ">GC\n";
	plotGC($svg,100,450,$fai,$chr,$start,$end,$g_step,$g_resolution);

	## mappability
	print ">mappability\n";
	plotMappability($svg,100,450,$MAP,$chr,$start,$end,$g_step,$g_resolution,$g_i);
	
	### gene
	print ">gene\n";
	my $hl=plotGene($svg,100,500,$GENE,$chr,$start,$end,$g_step,$g_resolution);
	
	## high light region
	print ">highlight region\n";
	plotHightlight($svg,100,200,$hl-200,\@posi,$chr,$start,$end,$g_step,$g_resolution);

	### output 
	print $out $svg->xmlify;
	close $out;
}
					
sub usage
{
	die `pod2text $0`;
}

sub plotPEM
{
	my ($svg,$plotX,$plotY,$in,$chr,$start,$end,$step,$resolution)=@_;
	my ($x1,$x2,$x3,$y);
	while(defined(my $line=<$in>)) 
	{
		chomp $line;
		#HWI-ST1328:55:D2A37ACXX:6:1308:16852:30958      5       141011678       5       141579202       567424
		my @colu=split(/\t/,$line);

		$x1=($colu[2]-$start)*$resolution/$step+$plotX;
		$x3=($colu[4]-$start)*$resolution/$step+$plotX;
		$y=$plotY-45;

		if(($colu[1] eq $chr) && ($colu[3] eq $chr) && ($colu[2]>=$start) && ($colu[2]<=$end) && ($colu[4]>=$start) && ($colu[4]<=$end) && ($colu[4]>=$colu[2])) 
		{
			my $leng=$colu[4]-$colu[2];
			if($leng>0) 
			{
				$svg->path('d',"M$x1,$plotY C$x1,$y $x3,$y $x3,$plotY",'style',{'fill',"none",'stroke',"orange"});
			}
		}else
		{
			if($colu[1] eq $chr && ($colu[2]>=$start) && ($colu[2]<=$end))
			{
				#left point
				$svg->line('x1',$x1,'y1',$plotY,'x2',$x1,'y2',$y-30,'stroke',"purple",'stroke-width',1);
			}
			if($colu[3] eq $chr && ($colu[4]>=$start) && ($colu[4]<=$end))
			{
				#right point
				$svg->line('x1',$x3,'y1',$plotY,'x2',$x3,'y2',$y-30,'stroke',"purple",'stroke-width',1);
			}
		}
	}
	my $txt_x=$plotX-50;
	my $txt_y=$plotY-15;
	$svg->text('x',0,'y',0,'stroke','black','fill','black','stroke-width',1,'-cdata',"PEM",'text-anchor','middle','transform',"translate($txt_x,$txt_y) rotate(-90)");
}
sub plotHightlight
{
	my ($svg,$plotX,$plotY,$plotHeight,$posi,$chr,$start,$end,$step,$resolution)=@_;
	my $x1=($posi->[0]-$start)*$resolution/$step+$plotX;
	my $x2=($posi->[1]-$start)*$resolution/$step+$plotX;
	#$svg->rect('x',$x1, 'y',$plotY,'width',$x2-$x1+1,'height',$plotHeight,'stroke','red','stroke-dasharrary','15 15','fill','pink','opacity',0.5);
	#$svg->text('x',($x1+$x2)/2,'y',$plotY-10,'stroke','black','fill','black','stroke-width',0.5,'-cdata',"$chr:$posi->[0]-$posi->[1]",'text-anchor','middle');
	if($start<=$posi->[0] && $posi->[0]<=$end && $start<=$posi->[1] && $posi->[1]<=$end)
	{
		$svg->rect('x',$x1, 'y',$plotY,'width',$x2-$x1+1,'height',$plotHeight,'stroke','red','stroke-dasharrary','15 15','fill','pink','opacity',0.5);
		$svg->text('x',($x1+$x2)/2,'y',$plotY-10,'stroke','black','fill','black','stroke-width',0.5,'-cdata',"$chr:$posi->[0]-$posi->[1]",'text-anchor','middle');
	}elsif($start<=$posi->[0] && $posi->[0]<=$end)
	{
		$svg->rect('x',$x1-1*$resolution*0.5, 'y',$plotY,'width',2*$resolution*0.5,'height',$plotHeight,'stroke','red','stroke-dasharrary','15 15','fill','pink','opacity',0.5);
		$svg->text('x',$x1,'y',$plotY-10,'stroke','black','fill','black','stroke-width',0.5,'-cdata',"$chr:$posi->[0]",'text-anchor','middle');
	}elsif($start<=$posi->[1] && $posi->[1]<=$end)
	{
		$svg->rect('x',$x2-1*$resolution*0.5, 'y',$plotY,'width',2*$resolution*0.5,'height',$plotHeight,'stroke','red','stroke-dasharrary','15 15','fill','pink','opacity',0.5);
		$svg->text('x',$x2,'y',$plotY-10,'stroke','black','fill','black','stroke-width',0.5,'-cdata',"$chr:$posi->[1]",'text-anchor','middle');
	}
}
sub plotAxis
{
	my ($svg,$plotX,$plotY,$plotHeight,$chr,$start,$end,$step,$resolution,$isRight)=@_;
	### X axis
	$svg->line('x1',$plotX,'y1',$plotY,'x2',($end-$start)*$resolution/$step+$plotX,'y2',$plotY,'stroke',"black",'stroke-width',1);
	### X axis's marker
	my $yu=$start%1000;
	my $yui=$start-$yu;
	my $nW=int(($end-$start)/1000);
	my $nInterval=int($nW/3);
	$yui+=$nInterval*1000;
	my ($x,$y,$mark);
	#print ">yui:$yui\t$chr\t$start\t$end\n";
	if($nInterval>1)
	{
		while($yui<$end) 
		{
			$x=($yui-$start)*$resolution/$step+$plotX;
			$svg->line('x1',$x,'y1',$plotY,'x2',$x,'y2',$plotY+5,'stroke',"black",'stroke-width',1);
			$mark=$yui/1000;
			$mark=$mark."K";
			$svg->text('x',$x,'y',$plotY+20,'stroke','black','fill','black','stroke-width',0.5,'text-anchor','middle','-cdata',$mark);
			$yui+=$nInterval*1000;
			#print "\tyui:$yui\tnInterval:$nInterval\n";
		}
	}
	### Y axis
	my $yInterval=int($plotHeight/2);
	if($isRight)
	{
		$svg->line('x1',$plotX+($end-$start)*$resolution/$step,'y1',$plotY,'x2',$plotX+($end-$start)*$resolution/$step,'y2',$plotY-$plotHeight,'stroke',"black",'stroke-width',1);
		foreach my $n (0..2) 
		{
			$y=$plotY-$yInterval*$n;
			$svg->line('x1',$plotX+($end-$start)*$resolution/$step,'y1',$y,'x2',$plotX+($end-$start)*$resolution/$step+5,'y2',$y,'stroke',"black",'stroke-width',1);
			$mark=$yInterval*$n/100;
			$svg->text('x',$plotX+($end-$start)*$resolution/$step+10,'y',$y,'stroke','black','fill','black','stroke-width',0.5,'text-anchor','left','-cdata',$mark);
		}
	}else
	{
		$svg->line('x1',$plotX,'y1',$plotY,'x2',$plotX,'y2',$plotY-$plotHeight,'stroke',"black",'stroke-width',1);
		foreach my $n (0..2) 
		{
			$y=$plotY-$yInterval*$n;
			$svg->line('x1',$plotX-5,'y1',$y,'x2',$plotX,'y2',$y,'stroke',"black",'stroke-width',1);
			$mark=$yInterval*$n;
			$svg->text('x',$plotX-30,'y',$y,'stroke','black','fill','black','stroke-width',0.5,'text-anchor','left','-cdata',$mark);
		}
	}
}
sub plotRD
{
	my ($svg,$plotX,$plotY,$in,$chr,$start,$end,$step,$resolution)=@_;
	## step: window size (default 100)
	my $n=0;
	my ($x1,$y1,$x2,$y2,$m1,$m2);
	my $ret=0;
	while(defined(my $line=<$in>)) 
	{
		chomp $line;
		##chrom  beg     end     i       RD_merge        RD_partition    RD_GC1  RD_GC2  RD_GC3  RD_raw
		#chr19   46805701        46805800        468057  31.58   30.65   30.56   30.33   30.65   34.00
		my @colu=split(/\t/,$line);
		my $cnv_chr=$chr;
		$colu[0]=~s/^chr//;
		#print STDERR "$cnv_chr\t$start\t$end\t$line\n";
		if($colu[1]>=$start && $colu[2]<=$end && $colu[0] eq $cnv_chr) 
		{
			$n++;
			## one pix per window
			#print STDERR "n:$n\tx1:$x1\ty1:$y1\n";
			if($n==1) {
				$x1=($colu[1]-$start)*$resolution/$step+$plotX;
				$y1=($plotY-$colu[4]);
				$m1=($plotY-$colu[9]);
			}else{
				$x2=($colu[1]-$start)*$resolution/$step+$plotX;
				$y2=($plotY-$colu[4]);
				$m2=($plotY-$colu[9]);
				#$svg->line('x1',$x1,'y1',$m1,'x2',$x2,'y2',$m2,'stroke',"red",'stroke-width',1,'opacity',0.8);
				#$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke',"green",'stroke-width',1);
				## step-like
				$svg->line('x1',$x1,'y1',$m1,'x2',$x1,'y2',$m2,'stroke',"red",'stroke-width',1,'opacity',0.8);
				$svg->line('x1',$x1,'y1',$m2,'x2',$x2,'y2',$m2,'stroke',"red",'stroke-width',1,'opacity',0.8);
				$svg->line('x1',$x1,'y1',$y1,'x2',$x1,'y2',$y2,'stroke',"green",'stroke-width',1,'opacity',0.8);
				$svg->line('x1',$x1,'y1',$y2,'x2',$x2,'y2',$y2,'stroke',"green",'stroke-width',1,'opacity',0.8);
				$x1=$x2;
				$y1=$y2;
				$m1=$m2;
			}
		}
	}
	print ">here\n";
	plotAxis($svg,$plotX,$plotY,100,$chr,$start,$end,$step,$resolution);
	print ">here2\n";
	my $txt_x=$plotX-50;
	my $txt_y=$plotY-50;
	$svg->text('x',0,'y',0,'stroke','black','fill','black','stroke-width',1,'-cdata',"Depth(X)",'text-anchor','middle','transform',"translate($txt_x,$txt_y) rotate(-90)");
	#$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','16','-cdata',"$sampleID",'text-anchor','middle','fill',"black",'transform',"translate($txt_x,$txt_y) rotate(-90)");
}

sub plotGC
{
	my ($svg,$plotX,$plotY,$fai,$chr,$start,$end,$step,$resolution)=@_;
	
	my $i=0;
	my $offset=$start%100;
	$offset=100-$offset;
	my $plotOffset=$offset*$resolution/$step;
	my $x0=$plotX+$plotOffset;
	my $y0=$plotY;
	my ($x1,$y1,$win,$gc);
	#my $_tStart=$start+$offset;
	my $target = $fai->fetch("$chr:$start-$end");
	while($step*$i+$offset<$end-$start+1)
	{
		# gc content
		#printf "offset:$offset\tstart:$start\tend:$end\tstep:$step\ti:$i\t%s\n",$step*$i+$offset;
		$win=substr($target,$step*$i+$offset,$step);
		$gc=0;
		while($win=~/G|C/gi) { $gc++; }
		$gc=int($gc*100/length($win));
		## one pix per window
		$x1=$plotX+$i*$resolution+$plotOffset;
		$y1=$plotY-$gc;
		if($i==1)
		{
			$svg->line('x1',$x0,'y1',$y1,'x2',$x1,'y2',$y1,'stroke',"blue",'stroke-width',1,'opacity',0.8);
		}
		if($i>=1)
		{
			$svg->line('x1',$x0,'y1',$y0,'x2',$x0,'y2',$y1,'stroke',"blue",'stroke-width',1,'opacity',0.8);
			$svg->line('x1',$x0,'y1',$y1,'x2',$x1,'y2',$y1,'stroke',"blue",'stroke-width',1,'opacity',0.8);
		}
		$x0=$x1;
		$y0=$y1;
		$i++;
	}
	plotAxis($svg,$plotX,$plotY,100,$chr,$start,$end,$step,$resolution);
	my $txt_x=$plotX-50;
	my $txt_y=$plotY-50;
	$svg->text('x',0,'y',0,'stroke','black','fill','black','stroke-width',1,'-cdata',"GC(%)",'text-anchor','middle','transform',"translate($txt_x,$txt_y) rotate(-90)");
}

sub plotMappability
{
	my ($svg,$plotX,$plotY,$in,$chr,$start,$end,$step,$resolution,$iLine)=@_;
	my $target="";
	my $_i=0;
	while(<$in>)
	{
		chomp;
		$_i++;
		if($_i == $iLine)
		{
			$target.=$_;
		}
	}
	my $offset=$start%100;
	$offset=100-$offset;
	my $plotOffset=$offset*$resolution/$step;
	my $x0=$plotX+$plotOffset;
	my $y0=$plotY;
	my ($x1,$y1,$win,$gc,$win_mappability);
	## step: window size (default 100)
	my $i=0;
	while($step*$i+$offset<$end-$start+1)
	{
		#print STDERR "i\t$i\tx0\t$x0\tx1\t$x1\n";
		$win=substr($target,$step*$i+$offset,$step);
		$win_mappability=0;
		my @_w=split //,$win;
		foreach (@_w) { $win_mappability+=$ENCODING{$_}; }
		$win_mappability=$win_mappability/(@_w);
		
		$x1=$plotX+$i*$resolution+$plotOffset;
		$y1=$plotY-$win_mappability*100;
		if($i==1)
		{
			$svg->line('x1',$x0,'y1',$y1,'x2',$x1,'y2',$y1,'stroke',"green",'stroke-width',1,'opacity',0.8);
		}
		if($i>=1)
		{
			$svg->line('x1',$x0,'y1',$y0,'x2',$x0,'y2',$y1,'stroke',"green",'stroke-width',1,'opacity',0.8);
			$svg->line('x1',$x0,'y1',$y1,'x2',$x1,'y2',$y1,'stroke',"green",'stroke-width',1,'opacity',0.8);
		}
		$x0=$x1;
		$y0=$y1;
		$i++;
	}
	plotAxis($svg,$plotX,$plotY,100,$chr,$start,$end,$step,$resolution,1);
	my $txt_x=$plotX+($end-$start)*$resolution/$step+50;
	my $txt_y=$plotY-50;
	$svg->text('x',0,'y',0,'stroke','black','fill','black','stroke-width',1,'-cdata',"Mappability",'text-anchor','middle','transform',"translate($txt_x,$txt_y) rotate(-90)");
}

sub plotGene
{
	my ($svg,$plotX,$plotY,$in,$chr,$start,$end,$step,$resolution)=@_;
	my $_iGene=-1;
	my $_yGene=0;
	my %hash=();
	my $n=0;
	my ($x1,$y1,$x2,$y2,$wid1,$wid2);
	my $retY=$plotY;
	while(defined(my $line=<$in>)) 
	{
		chomp $line;
		my  @colu=split(/\t/,$line);
		#970     NM_213590       chr13   +       50571142        50592603        50586076        50587300        2       50571142,50586070, 50571899,50592603,      0       TRIM13  cmpl    cmpl    -1,0,
		my $gene_chr=$chr;
		$colu[2]=~s/^chr//;
		if($gene_chr ne $colu[2]) { next; }
		# transcript overlap with the plot region
		if(($start<= $colu[4] && $colu[4]<=$end) || ($start<= $colu[5] && $colu[5]<=$end ) || ($colu[4]<= $start && $start<=$colu[5]) || ($colu[4]<= $end && $end<=$colu[5])) 
		{
			my $mark=$colu[12];
			if(!defined($hash{$mark}))
			{
				my @exon1=split(/,/,$colu[9]);
				my @exon2=split(/,/,$colu[10]);  
				$x2=$plotX+($colu[4]-$start)*$resolution/$step;
				$wid2=($colu[5]-$colu[4])*$resolution/$step;
				$_iGene++;
				$_yGene=$plotY+$_iGene*15;
				$retY=$_yGene+10;
				#{ print STDERR "$mark\t$chr\t$start\t$end\t$_yGene\t$_iGene\n"; }	
				$n++;
				$hash{$mark}=1;
				# whole transcript
				# strand
				$svg->line('x1',$x2>0?$x2:0,'y1',$_yGene,'x2',$x2+$wid2-1>0?$x2+$wid2-1:0,'y2',$_yGene,'stroke','blue');
				my $arrowCmd;
				my $strand=$colu[3];
				if($strand eq "-") { $arrowCmd="M 2,2 L 0,0 L 2,-2"; }
				else { $arrowCmd="M -2,2 L 0,0 L -2,-2"; }
				for(my $kk=$x2+10;$kk<$x2+$wid2-1;$kk+=10)
				{
					my $arrowY=$_yGene;
					$svg->path('d',$arrowCmd,'style',{'stroke',"blue",'fill','none'},'transform',"translate($kk,$arrowY) scale(1,1)");
				}
				$svg->text('x',$x2+$wid2+4,'y',$_yGene+6,'stroke','black','fill','black','stroke-width',1,'-cdata',$mark);
				#if($n%3==0)
				#{
				#	$svg->text('x',$x2,'y',530,'stroke','black','fill','black','stroke-width',1,'-cdata',$mark);
				#}elsif($n%3==1)
				#{
				#	$svg->text('x',$x2,'y',550,'stroke','black','fill','black','stroke-width',1,'-cdata',$mark);
				#}else
				#{
				#	$svg->text('x',$x2,'y',570,'stroke','black','fill','black','stroke-width',1,'-cdata',$mark);
				#}
				foreach $_ (0..$#exon1) 
				{
					$x1=$plotX+($exon1[$_]-$start)*$resolution/$step;
					$wid1=($exon2[$_]-$exon1[$_])*$resolution/$step;
					# exon
					$svg->rect('x',$x1?$x1:0, 'y',$_yGene-5,'width',$wid1,'height',10,'stroke','none','fill','blue');
				}
				#$x1=($colu[4]-$start)/$step+$plotX;
				#$wid1=($colu[6]-$colu[4])/$step;
				#$x2=($colu[7]-$start)/$step+$plotX;
				#$wid2=($colu[5]-$colu[7])/$step;
				# UTR
				#$svg->rect('x',$x1?$x1:0, 'y',$_yGene-7,'width',$wid1,'height',14,'stroke','none','fill','blue');
				#$svg->rect('x',$x2?$x2:0, 'y',$_yGene-7,'width',$wid2,'height',14,'stroke','none','fill','blue');
			}
	  }
	}
	$svg->text('x',$plotX-60,'y',$plotY-20,'stroke','black','fill','black','stroke-width',1,'-cdata',"Gene");		
	return $retY;
}

sub readEncoding
{
	#' '~[0-0]
	#'!'~[1-1]
	#'"'~[2-2]
	#'#'~[3-3]
	my ($pList,$infile)=@_;
	my $in;
	open $in,$infile or die "$!";
	while(<$in>)
	{
		chomp;
		/'(.)'~\[(\d+)\-(\d+)\]/;
		if($1 eq " ")
		{
			$pList->{$1}=-1;
		}else
		{
			$pList->{$1}=1/$2;
		}
	}
}
