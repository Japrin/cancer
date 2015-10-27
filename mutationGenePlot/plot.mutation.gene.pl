#!/usr/bin/perl
#============================================================================
# Name        		: plot.mutation.gene.pl
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

	perl plot.mutation.gene.pl [option] <out-prefix> <mutation (per gene)> <protein (tabix-ed gff)>
	
	-g	gene symbol [required]
	-h	display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use SVG;
use Data::Dumper;
use Tabix;
use List::Util qw(min max);
use URI::Escape;

my ($in,$out);
my ($opt_h,$opt_g);
GetOptions("h"	=>\$opt_h,"g=s"=>\$opt_g);
if(@ARGV<3 || $opt_h || !defined($opt_g)) { usage(); }
my $outfile=shift @ARGV;
my $mutfile=shift @ARGV;
my $featfile=shift @ARGV;

my %FT_Table=(
		#Molecule processing
		'Initiator methionine'=>{},
		'Signal peptide'=>{},
		'Transit peptide'=>{},
		'Propeptide'=>{},
		'Chain'=>{},
		'Peptide'=>{},
		## regions
		'Topological domain'=>{'section'=>'regions'},
		'Transmembrane'=>{'section'=>'regions'},
		'Intramembrane'=>{'section'=>'regions'},
		'Domain'=>{'section'=>'regions'},
		'Repeat'=>{'section'=>'regions'},
		'Calcium binding'=>{'section'=>'regions'},
		'Zinc finger'=>{'section'=>'regions'},
		'DNA binding'=>{'section'=>'regions'},
		'Nucleotide binding'=>{'section'=>'regions'},
		'Region'=>{'section'=>'regions'},
		'Coiled coil'=>{'section'=>'regions'},
		'Motif'=>{'section'=>'regions'},
		'Compositional bias'=>{'section'=>'regions'},
		## sites
		'Active site'=>{'section'=>'sites'},
		'Metal binding'=>{'section'=>'sites'},
		'Binding site'=>{'section'=>'sites'},
		'Site'=>{'section'=>'sites'},
		## amino acid 
		##		modifications
		'Non-standard residue'=>{'section'=>'aa_mod'},
		'Modified residue'=>{'section'=>'aa_mod'},
		'Lipidation'=>{'section'=>'aa_mod'},
		'Glycosylation'=>{'section'=>'aa_mod'},
		'Disulfide bond'=>{'section'=>'aa_mod'},
		'Cross-link'=>{'section'=>'aa_mod'},
		##		natural variatoins
		'Alternative sequence'=>{},
		'Natural variant'=>{},
		##		Experimental info
		'Mutagenesis'=>{},
		'Sequence uncertainty'=>{},
		'Sequence conflict'=>{},
		'Non-adjacent residues'=>{},
		'Non-terminal residue'=>{},
		## secondary structure
		'Helix'=>{'section'=>'2struc','color'=>'url(#myLinearGradient2)'},
		'Beta strand'=>{'section'=>'2struc','color'=>'url(#myLinearGradient3)'},
		'Turn'=>{'section'=>'2struc','color'=>'url(#myLinearGradient4)'},
		);

#star command
my $starCol="green";
my $starCmd="M -0.951,-0.309 L  0.951,-0.309 -0.588, 0.809 0.000,-1.000 0.588, 0.809 Z";
my $triangle1Cmd="M 0,0 L -5,-15 5,-15 Z";
my $triangle1Col="darkred";

### begin plot
my %mut=();
my %pnames=readMut(\%mut,$mutfile,$opt_g);
foreach my $pname (keys %pnames)
{
	my %gff=();
	printf STDERR "hehe\n";
	if(readGFF(\%gff,$featfile,$pname))
	{
		my $prot_len=$gff{'range'}->[1];
		
		my $margin=100;
		my $prot_h=30;
		my $prot_x=$margin;
		#margin + mutation region + title region
		my $prot_y=50+max(@{$mut{$pname}}*15,200)+100;
		my $wPerAA=2;
		my $width=max(2*$margin+$wPerAA*$prot_len,800);
		# prot_y + pro_h + annotation region
		my $height=max(($gff{'iRegions'}+3)*20,300)+$prot_h+$prot_y+$margin;
		#printf "$pname\twidth:$width\theight:$height\t$prot_y\t%d\n",max(($gff{'iRegions'}+3)*20,300);
		
		my $canvas = SVG->new('width',$width,'height',$height);
		$canvas->rect(x=>0,y=>0,width=>$width,height=>$height,fill=>"white");

		### drawing element
		#my $m=$canvas->marker('id','arrow','viewBox',"0 0 10 10",'refX',"0",'refY',"5",'markerUnits',"strokeWidth",'markerWidth',"8",'markerHeight',"6",'orient',"auto");
		#$m->path('d',"M 0 0 L 15 5 L 0 10 z");
		#my $m2=$canvas->marker('id','arrow2','viewBox',"0 0 10 10",'refX',"0",'refY',"5",'markerUnits',"strokeWidth",'markerWidth',"8",'markerHeight',"6",'orient',"auto",'transform',"rotate(180)");
		#$m2->path('d',"M 15 0 L 0 5 L 15 10 z");
		
		my $gfill1=$canvas->tag('radialGradient','id','myRadialGradient1','r','65%');
		$gfill1->tag('stop','offset','5%','stop-color','red');
		$gfill1->tag('stop','offset','95%','stop-color','darkred');
		my $gfill2=$canvas->tag('linearGradient','id','myLinearGradient2','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
		$gfill2->tag('stop','offset','5%','stop-color','dodgerblue');
		$gfill2->tag('stop','offset','30%','stop-color','lightskyblue');
		$gfill2->tag('stop','offset','50%','stop-color','lightskyblue');
		$gfill2->tag('stop','offset','95%','stop-color','dodgerblue');
		my $gfill3=$canvas->tag('linearGradient','id','myLinearGradient3','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
		$gfill3->tag('stop','offset','5%','stop-color','mediumseagreen');
		$gfill3->tag('stop','offset','30%','stop-color','palegreen');
		$gfill3->tag('stop','offset','50%','stop-color','palegreen');
		$gfill3->tag('stop','offset','95%','stop-color','mediumseagreen');
		my $gfill4=$canvas->tag('linearGradient','id','myLinearGradient4','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
		$gfill4->tag('stop','offset','5%','stop-color','hotpink');
		$gfill4->tag('stop','offset','35%','stop-color','lightpink');
		$gfill4->tag('stop','offset','50%','stop-color','lightpink');
		$gfill4->tag('stop','offset','95%','stop-color','hotpink');
		
		my $gfillBlack2Grey=$canvas->tag('linearGradient','id','myLinearGradientBlack2Grey','x1',0,'y1',0,'x2',1,'y2',0,'gradientTransform','rotate(90)');
		$gfillBlack2Grey->tag('stop','offset','5%','stop-color','grey');
		$gfillBlack2Grey->tag('stop','offset','35%','stop-color','lightgrey');
		$gfillBlack2Grey->tag('stop','offset','50%','stop-color','lightgrey');
		$gfillBlack2Grey->tag('stop','offset','95%','stop-color','grey');
		
		#protein name
		$canvas->text('x',100,'y',70,'font-family','TimesNewRoman','font-size','22','-cdata',"$opt_g($pname), $prot_len aa",'fill',"black");
		plotGene(\%gff,$canvas,$prot_x,$prot_y,$wPerAA,$prot_h);
		plotMutation($mut{$pname},$canvas,$prot_x,$prot_y,$wPerAA,$prot_h);
		plotLegend($canvas,$width-325,50,225);

		### output
		open $out,">","${outfile}.${pname}.svg" or die "$!";
		print $out $canvas->xmlify;
		close $out;
	}
}

############################################################################

sub plotMutation
{
	my ($pAry,$canvas,$x,$y,$wPerAA,$height)=@_;
	my $mutY=$y;
	my %_tList=();
	for(my $i=0;$i<@$pAry;$i++)
	{
		my ($aaPos,$aaChange)=@{$pAry->[$i]};
		my $_x=$x+$wPerAA*($aaPos-0.5);
		$_tList{$aaPos}++;
		my $_y=$mutY-($_tList{$aaPos}-1)*15;
		my $txt_x=$_x+6;
		my $txt_y=$_y-15;
		$canvas->path('d',$triangle1Cmd,'style',{'fill',$triangle1Col},'transform',"translate($_x,$_y) scale(1,1)");
		$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','12','-cdata',"$aaChange",'fill',"black",'transform',"translate($txt_x,$txt_y) rotate(0)");
		#printf "_x:$_x\t_y:$_y\taaPos:$aaPos\taaChange:$aaChange\n";
	}
}
sub plotLegend
{
	my ($canvas,$x,$y,$w)=@_;
	$canvas->rect('x',$x,'y',$y,'width',$w,'height',100,'fill','none','stroke','black','stroke-width',1.5);
	
	$canvas->rect('x',$x+5,'y',$y+5,'width',30,'height',30,'rx',5,'ry',5,'fill','url(#myLinearGradient2)');
	$canvas->text('x',$x+37,'y',$y+5+20,'font-family','TimesNewRoman','font-size','12','-cdata','Helix','fill','black');
	$canvas->rect('x',$x+67,'y',$y+5,'width',30,'height',30,'rx',5,'ry',5,'fill','url(#myLinearGradient3)');
	$canvas->text('x',$x+99,'y',$y+5+20,'font-family','TimesNewRoman','font-size','12','-cdata','Beta strand','fill','black');
	$canvas->rect('x',$x+166,'y',$y+5,'width',30,'height',30,'rx',5,'ry',5,'fill','url(#myLinearGradient4)');
	$canvas->text('x',$x+198,'y',$y+5+20,'font-family','TimesNewRoman','font-size','12','-cdata','Turn','fill','black');
	
	$canvas->rect('x',$x+5,'y',$y+40,'width',10,'height',30,'fill','darkorange','stroke','none');
	$canvas->text('x',$x+17,'y',$y+60,'font-family','TimesNewRoman','font-size','12','-cdata','sites','fill','black');
	$canvas->rect('x',$x+49,'y',$y+40,'width',10,'height',30,'fill','darkcyan','stroke','none');
	$canvas->text('x',$x+61,'y',$y+60,'font-family','TimesNewRoman','font-size','12','-cdata','aa_mod','fill','black');
	$canvas->rect('x',$x+102,'y',$y+40,'width',10,'height',30,'fill',"darkseagreen",'stroke','none');
	$canvas->text('x',$x+115,'y',$y+60,'font-family','TimesNewRoman','font-size','12','-cdata','region','fill','black');
	
	my $_x=$x+10;
	my $_y=$y+92;
	$canvas->path('d',$triangle1Cmd,'style',{'fill',$triangle1Col},'transform',"translate($_x,$_y) scale(1,1)");
	$canvas->text('x',$x+17,'y',$y+87,'font-family','TimesNewRoman','font-size','12','-cdata','sample(s)','fill','black');
	$_x=$x+77;
	$_y=$y+92;
	$canvas->path('d',$triangle1Cmd,'style',{'fill','lightgrey'},'transform',"translate($_x,$_y) scale(1,1)");
	$canvas->text('x',$x+84,'y',$y+87,'font-family','TimesNewRoman','font-size','12','-cdata','COSMIC','fill','black');
}
sub plotGene
{
	my ($pList,$canvas,$x,$y,$wPerAA,$height)=@_;
	my @range=@{$pList->{'range'}};
	#printf "range:\t%s\n",join("\t",@range);
	#printf "\$wPerAA: $wPerAA\t\$height: $height\n";
	#printf "$x<---------->%d\n",$wPerAA*$range[1];
	my $p=$x;
	#my $bgCol="lightcoral";
	my $bgCol="black";

	#$canvas->rect('x',$p,'y',$y,'width',$wPerAA*$range[1],'height',$height,'style',{'fill',$bgCol,'stroke','black'});
	$canvas->rect('x',$p,'y',$y,'width',$wPerAA*$range[1],'height',$height,'rx',5,'ry','5','style',{'fill','url(#myLinearGradientBlack2Grey)','stroke','black'});
	my $txtN_x=$p-15;
	my $txtN_y=$y+$height/2+4;
	my $txtC_x=$p+$wPerAA*$range[1]+4;
	my $txtC_y=$txtN_y;
	$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','12','-cdata',"N",'fill',"black",'transform',"translate($txtN_x,$txtN_y) rotate(0)");
	$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','12','-cdata',"C",'fill',"black",'transform',"translate($txtC_x,$txtC_y) rotate(0)");

	my $jj=3;
	foreach my $pFT(@{$pList->{'feature'}})
	{
		# feature, beg, end, note
		# feat	beg	end	note
		my ($_feat,$_beg,$_end,$_note)=@$pFT;
		my $_x=$p+$wPerAA*($pFT->[1]-1);
		my $_w=$wPerAA*($pFT->[2]-$pFT->[1]+1);
		my $_h=10;
		my $_x2=$_x+$_w;
		my $_xMid=$_x+$_w/2;
		if($FT_Table{$pFT->[0]}->{'section'} eq '2struc')
		{
			#printf "2struc:\t%s\n",join("\t",@$pFT);
			my $_col=$FT_Table{$pFT->[0]}->{'color'};
			$canvas->rect('x',$_x,'y',$y,'width',$_w,'height',$height,'rx',3,'ry','3','style',{'fill',$_col,'opacity',0.8});
		}
		# other: regions, sites, aa_mod
		## P04637  UniProtKB       Helix   36      38      .       .       .       .
		
		my $_col="";
		my $_y=0;	
		my $_txt_y=0;
		my $_txt_x=0;
		my $_txt="";
		if($FT_Table{$pFT->[0]}->{'section'} eq 'regions') 
		{ 
			$_col="darkseagreen";
			$_y=$y+10+$height+$jj*20;
			$_txt_y=$_y+8-10;
			$_txt_x=$_x+$_w/2;
			$_txt="$pFT->[0]: $_note";
			$canvas->rect('x',$_x,'y',$_y,'width',$_w,'height',$_h,'style',{'fill',$_col,'stroke','none','opacity',0.80});
			$canvas->text('x',0,'y',0,'font-family','TimesNewRoman','font-size','10','-cdata',"$_txt",'fill',"black","text-anchor","middle",'transform',"translate($_txt_x,$_txt_y) rotate(0)");
			$jj++; 
		}elsif($FT_Table{$pFT->[0]}->{'section'} eq 'sites')
		{
			$_col="darkorange"; 
			$_y=$y+5+$height+1*20;
			$canvas->rect('x',$_x,'y',$_y,'width',$_w,'height',$_h,'style',{'fill',$_col,'stroke','none','opacity',0.80});
		}elsif($FT_Table{$pFT->[0]}->{'section'} eq 'aa_mod')
		{
			$_col="darkcyan"; 
			$_y=$y+5+$height+2*20;
			if($pFT->[0] eq "Disulfide bond")
			{
				$canvas->rect('x',$_x,'y',$_y,'width',$wPerAA,'height',$_h,'style',{'fill',$_col,'stroke','none','opacity',0.80});
				$canvas->rect('x',$_x2,'y',$_y,'width',$wPerAA,'height',$_h,'style',{'fill',$_col,'stroke','none','opacity',0.80});
			}else
			{
				$canvas->rect('x',$_x,'y',$_y,'width',$_w,'height',$_h,'style',{'fill',$_col,'stroke','none','opacity',0.80});
			}
		}
		
		#printf "$FT_Table{$pFT->[0]}->{'section'}:\t%s\t_y:$_y\n",join("\t",@$pFT);
	}
}
sub readMut
{
	my ($pList,$infile,$gname)=@_;
	my $in;
	open $in,$infile or die "$!";
	my %pnames=();
	while(<$in>)
	{
		chomp;
		my @F=split /\t/;
		if($F[3] ne $gname) { next; }
		#NM_000546       chr17   7577538 TP53    exon7   c.G743A p.R248Q NP_000537       TP53    7157    K7PPA8
		$pnames{$F[-1]}=1;
		if(!exists($pList->{$F[-1]})) { $pList->{$F[-1]}=[]; }
		my ($pos)=$F[6]=~/(\d+)/;
		$F[6]=~s/p\.//;
		push @{$pList->{$F[-1]}},[$pos,$F[6]];
	}
	return %pnames;
}
sub readGFF
{
	my ($pList,$infile,$pname)=@_;
	my @feat=();
	my @range=();

	my $t = Tabix->new(-data=>"$infile");
	my $iter = $t->query($pname);
	my $iRegion=0;
	while($_=$t->read($iter))
	{
		my @F=split /\t/,$_;
		if($F[2] eq "RANGE")
		{
			@range=@F[3,4];
			next;
		}
		my $section=$FT_Table{$F[2]}->{'section'};
		if(!$section) { next; }
		if(!($section eq 'regions' || $section eq 'sites' || $section eq 'aa_mod' || $section eq '2struc')) { next; }
		## P04637  UniProtKB       Helix   36      38      .       .       .       .
		my ($note)=$F[8]=~/Note=(.+)/;
		$note  = uri_unescape($note);
		if(!$note) { $note=""; }
		my @p=(@F[2,3,4],$note);
		push @feat,\@p;
		if($section eq 'regions') { $iRegion++; }
	}
	if(@feat<1) { return 0; }
	@feat=sort {
				if($FT_Table{$a->[0]}->{'section'} eq "regions" && $FT_Table{$b->[0]}->{'section'} eq "regions") { $a->[1] <=> $b->[1]; }
				elsif($FT_Table{$a->[0]}->{'section'} eq "regions") { 1; }
				elsif($FT_Table{$b->[0]}->{'section'} eq "regions") { -1; }
				else {
					if(($FT_Table{$a->[0]}->{'section'} cmp $FT_Table{$b->[0]}->{'section'}) !=0 ) { $FT_Table{$a->[0]}->{'section'} cmp $FT_Table{$b->[0]}->{'section'}; }
					else { $a->[1] <=> $b->[1]; }
				}
			} @feat;
	$pList->{'feature'}=\@feat;
	$pList->{'range'}=\@range;
	$pList->{'iRegions'}=$iRegion;
	return 1;
}

sub usage
{
	die `pod2text $0`;
}
