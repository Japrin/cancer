#!/usr/bin/perl



use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;
use Term::ANSIColor;

our $REVISION = '$Revision: 5b1827d17ee9dcde01fd84cbd230edc3b98f9595 $';
our $DATE =	'$Date: 2015-03-22 15:29:51 -0700 (Sun, 22 Mar 2015) $';  
## modify here:
our $AUTHOR =	'$Author: Kai Wang <kai@openbioinformatics.org>; Liangtao Zheng <zhengliangtao@pku.edu.cn> $';
######our $AUTHOR =	'$Author: Kai Wang <kai@openbioinformatics.org> $';

our ($verbose, $help, $man, $mrnaseq,$opt_independent,$opt_protein,$opt_k,$opt_peptide,$opt_sample,$opt_vcf);

our ($evffile, $genefile, $fastafile, $out_protein, $out_peptide);

our %codon1 = (TTT=>"F", TTC=>"F", TCT=>"S", TCC=>"S", TAT=>"Y", TAC=>"Y", TGT=>"C", TGC=>"C", TTA=>"L", TCA=>"S", TAA=>"*", TGA=>"*", TTG=>"L", TCG=>"S", TAG=>"*", TGG=>"W", CTT=>"L", CTC=>"L", CCT=>"P", CCC=>"P", CAT=>"H", CAC=>"H", CGT=>"R", CGC=>"R", CTA=>"L", CTG=>"L", CCA=>"P", CCG=>"P", CAA=>"Q", CAG=>"Q", CGA=>"R", CGG=>"R", ATT=>"I", ATC=>"I", ACT=>"T", ACC=>"T", AAT=>"N", AAC=>"N", AGT=>"S", AGC=>"S", ATA=>"I", ACA=>"T", AAA=>"K", AGA=>"R", ATG=>"M", ACG=>"T", AAG=>"K", AGG=>"R", GTT=>"V", GTC=>"V", GCT=>"A", GCC=>"A", GAT=>"D", GAC=>"D", GGT=>"G", GGC=>"G", GTA=>"V", GTG=>"V", GCA=>"A", GCG=>"A", GAA=>"E", GAG=>"E", GGA=>"G", GGG=>"G");

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 
	'mrnaseq'=>\$mrnaseq, 
	'independent|i'=>\$opt_independent,
	'protein|p=s'=>\$opt_protein,
	'kmer|k=i'=>\$opt_k,
	'peptide|e=s'=>\$opt_peptide,
	'sample|s=s'=>\$opt_sample,
	'vcf'=>\$opt_vcf
) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 3 or pod2usage ("Syntax error");

($evffile, $genefile, $fastafile) = @ARGV;

if($opt_vcf && $evffile =~ /\.gz$/) { open EVF, "bgzip -cd $evffile | " or die "Error: cannot read from evffile $evffile: $!\n"; }
else { open (EVF, $evffile) or die "Error: cannot read from evffile $evffile: $!\n"; }
open (GENE, $genefile) or die "Error: cannot read from genefile $genefile: $!\n";
open (FASTA, $fastafile) or die "Error: cannot read from fastafile $fastafile: $!\n";

if(defined($opt_protein)) { open $out_protein,">",$opt_protein or die "$!"; }
if(!defined($opt_k)) { $opt_k=9; }
if(defined($opt_peptide)) { open $out_peptide,">",$opt_peptide or die "$!"; }
if(!defined($opt_sample)) { $opt_sample=""; }

### mutation type definition
use constant {
	SINGLEBASE_SUBSTITUTION => 0,
	SINGLEBASE_DELETION => 1,
	MULTIBASE_DELETION => 2,
	BLOCK_SUBSTITUTION => 3,
	INSERTION => 4,
	};
	

## @ read in variant list;
## @ important info: transcript, start, end, obs
## modify here:
my (%var_trans, %need_trans);
######my (@queue, %need_trans);
while (<EVF>) {
	if(/^#/) { next; }
	s/[\r\n]+$//;
	####m/^line\d+/ or die "Error: invalid record found in exonic_variant_function file $evffile (line number expected): <$_>\n";
	my @field = split (/\t/, $_);
	my ($evfInfo,$lineInfo,$evType);
	if($opt_vcf)
	{
		($evfInfo)=/AAChange=(.+?);/;
		$lineInfo="$field[0]:$field[1]";
		($evType)=/ExonicFunc=(.+?);/;
		my ($vType)=/Func=(.+?);/;
		if($vType ne "exonic" && $vType ne "exonic,splicing") { next; }
	}else
	{
		$evfInfo=$field[2];
		$lineInfo=$field[0];
		$evType=$field[1];
	}
	$evfInfo =~ m/^\w+:(\w+):wholegene/ and next;
	$evfInfo =~ m/unknown/i and next;
	
	## line3   nonsynonymous SNV       DDX11L1:uc010nxq.1:exon3:c.A187G:p.T63A,        1       13424   13424   A       G
	#ensembl gene name may contain - and ., where as UCSC transcript name may contain '.'.
	#$field[2] =~ m/^[\w\-\.]+?:([\w\.]+?):exon\d+:c.(\w+)/ or warn "Error: invalid record found in exonic_variant_function file (exonic format error): <$_>\n" and next;
	#$field[2] =~ m/^([\w\-\.]+?):([\w\.]+?):exon\d+:c.(\w+)(?::p.(\w+))*,/ or warn "Error: invalid record found in exonic_variant_function file (exonic format error): <$_>\n" and next;
	$evfInfo =~ m/^(.+?):([\w\.]+?):exon\d+:c.(\w+)(?::p.(.+))*\b/ or warn "Error: invalid record found in exonic_variant_function file (exonic format error): <$_>\t$evfInfo\n" and next;
	my ($gene, $transcript, $cchange, $aaChange) = ($1, $2, $3, $4);
	if($aaChange) { $aaChange=(split /,/,$aaChange)[0]; }
	my ($gChr,$gBeg,$gEnd,$gRef,$gAlt)=@field[3..7];
	my $gMut="$gChr:$gBeg-$gEnd:$gRef:$gAlt";

	if(!defined($aaChange)) { $aaChange=""; }
	my ($start, $end, $ref, $obs);
	
	if ($cchange =~ m/^(\w)(\d+)(\w)$/) {
		($start, $end, $ref, $obs) = ($2, $2, $1, $3);
	} elsif ($cchange =~ m/^(\d+)_(\d+)delins(\w+)/) {	#block substitution
		($start, $end, $obs) = ($1, $2, $3);
	} elsif ($cchange =~ m/^(\d+)del\w+/) {		#single base deletion
		($start, $end, $obs) = ($1, $1, '');
	} elsif ($cchange =~ m/^(\d+)_(\d+)del(\w*)/) {	#multi-base deletion
		($start, $end, $obs) = ($1, $2, '');
	} elsif ($cchange =~ m/^(\d+)_(\d+)ins(\w+)/) {	#insertion
		($start, $end, $ref, $obs) = ($1, $1, 0, $3);		#if end is equal to start, this is an insertion
	} elsif ($cchange =~ m/^(\d+)dup(\w+)/) {	#insertion
		($start, $end, $ref, $obs) = ($1, $1, 0, $2);
	} elsif ($cchange =~ m/^(\d+)_(\d+)(\w+)/) {	#non-frameshift substitution
		($start, $end, $obs) = ($1, $2, $3);
	} else {
		die "Error: invalid coding change format: <$cchange> within <$_>\n";
	}
	## modify here:
	push @{$var_trans{$transcript}},[$lineInfo, $transcript, $start, $end, $ref, $obs, $cchange, $aaChange, $gene, $gMut, $evType];
	######push @queue, [$field[0], $transcript, $start, $end, $ref, $obs, $cchange];
	$need_trans{$transcript}++;
}

## @ read in gene model file (genePred format);
## @ return mrnastart and mrnaend coordinate(in transcript coordinate system) for each transcript;
my (%mrnastart, %mrnaend);
while (<GENE>) {
	s/[\r\n]+$//;
	my @field = split (/\t/, $_);
	@field >= 11 or die "Error: invalid record found in gene file (>=11 fields expected): <$_>\n";
	$field[0] =~ m/^\d+$/ and shift @field;		#refGene and ensGene has bin as the first column

	my ($name, $strand, $txstart, $txend, $cdsstart, $cdsend, $exonstart, $exonend) = @field[0, 2, 3, 4, 5, 6, 8, 9];
	$need_trans{$name} or next;
	
	my ($mrnastart, $mrnaend);
	
	#next we need to make sure that there is no intron between transcription start and translation start (this is rare but it happens when cdsstart is not in the first exon)
	my @exonstart = split (/,/, $exonstart);
	my @exonend = split (/,/, $exonend);
	
	$txstart++;
	$cdsstart++;
	@exonstart = map {$_+1} @exonstart;

	if ($strand eq '+') {
		#<---->-----<--->----<------>----<----->-----<--->
		#             **********
		my $intron = 0;
		for my $i (0 .. @exonstart-1) {
			$i and $intron += ($exonstart[$i]-$exonend[$i-1]-1);
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				$mrnastart = $cdsstart-$txstart+1-$intron;
			}
			if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
				$mrnaend = $cdsend-$txstart+1-$intron;
			}
			
		}
	} elsif ($strand eq '-') {
		#<---->-----<--->----<------>----<----->-----<--->
		#             **********
		my $intron = 0;
		for (my $i=@exonstart-1; $i>=0; $i--) {
			$i<@exonstart-1 and $intron += ($exonstart[$i+1]-$exonend[$i]-1);
			if ($cdsend >= $exonstart[$i] and $cdsend <= $exonend[$i]) {
				$mrnastart = $txend-$cdsend+1-$intron;
			}
			if ($cdsstart >= $exonstart[$i] and $cdsstart <= $exonend[$i]) {
				$mrnaend = $txend-$cdsstart+1-$intron;
			}
			
		}
	}
			

	$mrnastart{$name} = $mrnastart;
	$mrnaend{$name} = $mrnaend;
}


## @ read in transcriptome file (FA format);
## @ 
my (%mrnaseq);
my ($curname, $curseq);

while (<FASTA>) {
	s/[\r\n]+$//;
	if (m/^>([\w\.]+)/) {
		if ($curseq) {
			$mrnaseq{$curname} = $curseq;
		}
		$curname = $1;
		$curseq = '';
	} else {
		$curseq .= $_;
	}
	$curseq and $mrnaseq{$curname} = $curseq;	#process the last sequence
}

#process each element in each mutated transcript
foreach my $transcript (keys %var_trans)
{
	my @_ary=@{$var_trans{$transcript}};
	my $gene=$var_trans{$transcript}->[0]->[8];
	
	## check
	if (not defined $mrnaseq{$transcript}) {
		print STDERR "WARNING: cannot find mRNA sequence for $transcript in the fastafile $fastafile\n";
		next;
	}
	if (not defined $mrnastart{$transcript}) {
		print STDERR "WARNING: cannot find annotation for $transcript in the genefile $genefile or cannot infer the transcription start site\n";
		next;
	}
	if (not defined $mrnaend{$transcript}) {
		print STDERR "WARNING: cannot find annotation for $transcript in the genefile $genefile or cannot infer the transcription end site\n";
		next;
	}
	
	## @dna_tx for the transcript's dna sequence in gene model file
	my $dna_tx = substr ($mrnaseq{$transcript}, $mrnastart{$transcript}-1, $mrnaend{$transcript}-$mrnastart{$transcript}+1);
	my @dna_tx = split (//, $dna_tx);

	my ($protein1, $protein2);
	my $warning = '';
	#$verbose and print STDERR "NOTICE: wild-type DNA sequence is $dna\n";
	$protein1 = translateDNA ($dna_tx);
	my $protein1_hasStop=1;
	if($protein1!~/\*/) { $protein1_hasStop=0; }
	$protein1=~s/\*$//;
	$protein1 =~ m/\*\w/ and $warning = '(WARNING: Potential FASTA sequence error!!!)';
	if($opt_independent)
	{
		#mutateSequenceArray($start,$end,$ref,$obs,$pAryDNA,$pAryDNAW,$pAryDNAWDis);
		for my $i (0 .. @_ary-1)
		{
			## $dna is what we will operate
			my @dna = @dna_tx;
			## @dna_w for wildtype dna; @dna_w_dis for wildtype dna's display only;
			my @dna_w=@dna_tx;
			my @dna_w_dis=@dna_tx;

			my ($line, $_tx, $start, $end, $ref, $obs, $cchange, $aaChange, $_gene, $gMut, $evType) = @{$_ary[$i]};
			if($evType=~/^synonymous/) { next; }

			if ($end > length ($mrnaseq{$transcript})) 
			{
				print STDERR "ERROR: transcript end ($mrnaend{$transcript}) for $transcript is longer than transcript length ${\(length ($mrnaseq{$transcript}))}, skipping this transcript\n";
				next;
			}
			if ($end > @dna) {
				print STDERR "ERROR in $line: end position of variant ($end) in $transcript is longer than coding portion length ${\(scalar @dna)}, skipping this transcript\n";
				next;
			}
			my ($aastart,$aaend1,$aaend2,$mutType)=mutateSequenceArray($start,$end,$ref,$obs,\@dna,\@dna_w,\@dna_w_dis);
			$protein2 = translateDNA (join('',@dna));
			### start codon loss
			if(substr($protein2,0,1) ne "M") { next; }

			my $isFrameShift=0;
			my $isStoploss=0;
			my $isStopgain=0;
			if($protein2!~/\*/)
			{
				$isStoploss=1;
			}else
			{
				$protein2=~s/\*$//;
				if($protein2=~/\*/) { $isStopgain=1; }
				$protein2 =~ s/\*.*//;
			}

			my $protein_id=sprintf "$gene.$transcript.Mutant%04d.$opt_sample",$i+1;
			my $protein_id_comment=sprintf "cDNAChange: [%s], aaChange: [%s], varSite alleles: [%s], mutType: [%s], genomeMut: [%s]", $cchange, $aaChange, ($dna[$start-1] eq ""?"-":$dna[$start-1]), formatMutType($mutType), $gMut;
			
			if($opt_protein)
			{	
				#print $out_protein ">$transcript.WILDTYPE\n";
				#print $out_protein "$protein1\n";
				printf $out_protein ">$protein_id $protein_id_comment\n";
				print $out_protein "$protein2\n";
			}
			if($verbose)
			{
				print STDERR colored(">$transcript.WILDTYPE\n","bold blue");
				my ($vis_cDNA,$vis_prot)=verbose_output_sequence(\@dna_w_dis,[$start-1],1);
				print STDERR "$vis_cDNA\n";
				print STDERR "$vis_prot\n";
				printf STDERR colored(">$protein_id $protein_id_comment\n","bold blue");
				($vis_cDNA,$vis_prot)=verbose_output_sequence(\@dna,[$start-1],1);
				print STDERR "$vis_cDNA\n";
				print STDERR "$vis_prot\n";
				print STDERR "\n";
			}
			if($out_peptide)
			{
				my @peptideListWildtype=makeKMerList($protein1,$opt_k,1,length($protein1),{});
				my %peptideListWildtype=();
				foreach (@peptideListWildtype) { $peptideListWildtype{$_->[0]}=1; }
				my @peptideList=();
				#if($isFrameShift)
				#if($isStopgain || $isStoploss)
				if($isStopgain || ($isStoploss && $protein1_hasStop))
				{
					@peptideList=makeKMerList($protein2,$opt_k,$aastart,length($protein2),\%peptideListWildtype);

				}else
				{
					@peptideList=makeKMerList($protein2,$opt_k,$aastart,$aaend2,\%peptideListWildtype,$mutType==MULTIBASE_DELETION);
				}

				for(my $_j=0;$_j<@peptideList;$_j++)
				{
					printf $out_peptide ">${protein_id}_Pep%04d %s/%s %s\n",$_j+1,$peptideList[$_j]->[1],$opt_k,$protein_id_comment;
					print $out_peptide "$peptideList[$_j]->[0]\n";
				}
			}
		}
	
	}else
	{
		## $dna is what we will operate
		my @dna = @dna_tx;
		## @dna_w for wildtype dna; @dna_w_dis for wildtype dna's display only;
		my @dna_w=@dna_tx;
		my @dna_w_dis=@dna_tx;

		my %varPos=();
		my @mutPos=();
		## use all variant in this transcript
		for my $i (0 .. @_ary-1)
		{
			my ($line, $_tx, $start, $end, $ref, $obs, $cchange, $aaChange, $_gene, $gMut, $evType) = @{$_ary[$i]};
			if ($end > length ($mrnaseq{$transcript})) 
			{
				print STDERR "ERROR: transcript end ($mrnaend{$transcript}) for $transcript is longer than transcript length ${\(length ($mrnaseq{$transcript}))}, skipping this transcript\n";
				next;
			}
			if ($end > @dna) {
				print STDERR "ERROR in $line: end position of variant ($end) in $transcript is longer than coding portion length ${\(scalar @dna)}, skipping this transcript\n";
				next;
			}

			my ($aastart,$aaend1,$aaend2,$mutType)=mutateSequenceArray($start,$end,$ref,$obs,\@dna,\@dna_w,\@dna_w_dis);
			$varPos{$start-1}=[$start, $end, $ref, $obs, $cchange, $aaChange, $mutType, $_gene, $gMut];

		}
		
		my @varPos=sort { $a<=>$b } keys %varPos;
		my @refAlleles=@dna_w[@varPos];
		my @varAlleles=@dna[@varPos];
		my @r=getAllPossibleCombination(\@refAlleles,\@varAlleles,scalar @refAlleles);
	
		if($verbose)
		{
			print STDERR "##\@varPos (0-based): ".join(" ",@varPos)."\n";
			print STDERR "##\@refAlleles: ".join(" ",map { ($_ eq '')?'-':$_ } @refAlleles)."\n";
			print STDERR "##\@varAlleles: ".join(" ",map { ($_ eq '')?'-':$_ } @varAlleles)."\n";
			print STDERR ">$transcript.WILDTYPE\n";
			my ($vis_cDNA,$vis_prot)=verbose_output_sequence(\@dna_w_dis,\@varPos,0);
			print STDERR "$vis_cDNA\n";
			#print STDERR "$vis_prot\n";
			print STDERR ">$transcript.Mutant (all possible combination: ".scalar @r.")\n";
			($vis_cDNA,$vis_prot)=verbose_output_sequence(\@dna,\@varPos,0);
			print STDERR "$vis_cDNA\n";

		}
		my @tx_out=();	
		my $j=0;
		foreach (@r)
		{
			my @_t=@dna;
			@_t[@varPos]=@$_;
			my $_protein_t = translateDNA (join('',@_t));
			
			my $isFrameShift=0;
			my $isStoploss=0;
			my $isStopgain=0;
			if($_protein_t!~/\*/)
			{
				$isStoploss=1;
			}else
			{
				$_protein_t=~s/\*$//;
				if($_protein_t=~/\*/) { $isStopgain=1; }
				$_protein_t =~ s/\*.+//;
			}
			
			### start codon loss or not
			if(substr($_protein_t,0,1) eq "M")
			{
				$j++;

				my $protein_id=sprintf "$gene.$transcript.Mutant%04d.$opt_sample",$j;
				my $protein_id_comment=sprintf "cDNAChange: [%s], aaChange: [%s], varSite alleles: [%s], mutType: [%s], genomeMut: [%s]", 
					join(" ",map { $varPos{$_}->[4] } @varPos),
					join(" ",map { $varPos{$_}->[5] } @varPos), 
					join(" ", map { ($_ eq '')?'-':$_ } @$_),
					join(" ",map { formatMutType($varPos{$_}->[6]) } @varPos),
					join(" ", map { $varPos{$_}->[8] } @varPos);
				push @tx_out,[$protein_id,$protein_id_comment,$_protein_t,$isStopgain,$isStoploss];

				if($verbose)
				{
					print STDERR colored(">$protein_id $protein_id_comment\n","bold blue");
					my ($vis_cDNA,$vis_prot)=verbose_output_sequence(\@_t,\@varPos,1);
					print STDERR "$vis_cDNA\n";
					print STDERR "$vis_prot\n";
					print STDERR "\n";
				}
			}
		}

		if($opt_protein)
		{
			foreach (@tx_out)
			{	
				#print $out_protein ">$transcript.WILDTYPE\n";
				#print $out_protein "$protein1\n";
				print $out_protein ">$_->[0] $_->[1]\n";
				print $out_protein "$_->[2]\n";
			}
		}
		if($out_peptide)
		{
			my @peptideListWildtype=makeKMerList($protein1,$opt_k,1,length($protein1),{});
			my %peptideListWildtype=();
			foreach (@peptideListWildtype) { $peptideListWildtype{$_->[0]}=1; }
			
			foreach my $_prot (@tx_out)
			{
				my @peptideList=();
				####my $isFrameShift=$_prot->[3];
				my $protein_id=$_prot->[0];
				my $protein_id_comment=$_prot->[1];
				@peptideList=makeKMerList($_prot->[2],$opt_k,1,scalar @dna,\%peptideListWildtype);

				for(my $_j=0;$_j<@peptideList;$_j++)
				{
					printf $out_peptide ">${protein_id}_Pep%04d %s/%s %s\n",$_j+1,$peptideList[$_j]->[1],$opt_k,$protein_id_comment;
					print $out_peptide "$peptideList[$_j]->[0]\n";
				}

			}

		}
	}
	
}

## @ given an array of sequence, and positins to be highlighted, return colored string
## @ pSeq, reference of the array of sequence
## @ pPos, postions to be highlighted, should be sorted increasingly
sub verbose_output_sequence
{
	my ($pSeq,$pPos,$toTranslate)=@_;
	my %h=();
	map { $h{$_}=1 } @$pPos;
	my @s=map { ($_ eq '')?'-':$_ } @$pSeq;
	###
	if(!$toTranslate)
	{
		#my @test=@s;
		for(my $i=0;$i<@$pSeq;$i++)
		{
			if($h{$i})
			{
				$s[$i]=colored($s[$i],'bold red');
			}else
			{
				$s[$i]=colored($s[$i],'green');
			}
		}
		return join(" ",@s);
	}else
	{
		my %hh=();
		### calculate positions in @s which corresponse to @$pPos
		my $k=0;
		## @ss for visualization, so with escaped strings;
		## @tt for calculation, so without escaped strings;
		my @ss=split //,join(" ",@s);
		my @tt=@ss;
		for(my $j=0;$j<@s;$j++)
		{
			# from k to k+length($s[$j])-1;
			if($h{$j})
			{
				for(my $_m=$k;$_m<=$k+length($s[$j])-1;$_m++)
				{
					$hh{$_m}=1;
				}
			}
			## s[j] length
			$k=$k+length($s[$j])-1;
			## a space
			$k+=2;
		}

		my @aa=();
		my $_phase=0;
		for(my $j=0;$j<@tt;$j++)
		{
			if($tt[$j] ne " " && $tt[$j] ne '-') 
			{ 
				$_phase++; 
			}
			if($_phase !=0 )
			{
				if($hh{$j})
				{
					$ss[$j]=colored($tt[$j],'underline bold red');
				}else
				{
					$ss[$j]=colored($tt[$j],'underline green');
				}
			}
			if($_phase ==1 && $tt[$j] ne " " && $tt[$j] ne '-')
			{
				## try translate the codon
				my $c=$tt[$j];
				my $l=$j+1;
				## j .. l determin on codon c
				for(;$l<@tt;$l++)
				{
					if($tt[$l] ne " " && $tt[$l] ne "-")
					{
						$c=$c.$tt[$l];
						#print STDERR "$j\t$l\t$c\n";
					}
					if(length($c)==3)
					{
						#print STDERR "$codon1{$c}\n";
						last;
					}
				}
				if(length($c)==3)
				{
					push @aa,$codon1{$c}.(" "x($l-$j));
					for(my $ll=$l+1;$ll<@tt;$ll++)
					{
						if($tt[$ll] eq " ")
						{
							push @aa," ";
						}elsif($tt[$ll] eq "-")
						{
							push @aa,"-";
						}else
						{
							last;
						}
					}
				}
			}
			
			if($_phase==3) { $_phase=0; }
		}
		for(my $j=0;$j<@aa;$j++)
		{
			$aa[$j]=colored($aa[$j],'green');
		}
		return (join("",@ss),join("",@aa));
	}
	
}

## @ give string $s, make "overlap"-mers constructed by concacate all kmers overlap positions from i to j, highlight position i and j (relative position in kmer )
## @ $k, $i, $j, is 1-based
## @ $isMultiBaseNonFrameShiftDel
## @ makeOverlapMer($s,$k,$i,$j)
sub makeOverlapMer
{
	my ($s,$k,$i,$j,$bgPepList,$isMultiBaseNonFrameShiftDel)=@_;
	if(!defined($isMultiBaseNonFrameShiftDel)) { $isMultiBaseNonFrameShiftDel=0; }
	my @s=split //,$s;
	my @ret=();
	## --------------------XXXX-----------------
	##                     i  j
	##             =========
	##             t
	##              =========
	##               =========
	##                        =========
	for(my $t=$i-($k-1);$t<=($j-$isMultiBaseNonFrameShiftDel) && $t<=@s-($k-1);$t++)
	{
		if($t<1) { next; }
		my $m="";
		if($i==$j) { $m=$i-$t+1; }
		else 
		{
		  	my $p1=$i-$t+1;
			my $p2=$j-$t+1;
			if($p2>$k) 
			{ 
				$p2=$k;
			}
			if($p1<1)
			{
				$p1=1;
			}
			if($p1==$p2) { $m="$p1"; }
			else { $m="$p1..$p2"; }
		}
		my $pep=join('',@s[$t-1 .. $t+$k-2]);
		if(!$bgPepList->{$pep})
		{
			push @ret,[$pep,$m];
		}
		
	}
	return @ret;
}

## @ give string $s, make all kmers overlap positions from i to j, highlight position i and j (relative position in kmer )
## @ $k, $i, $j, is 1-based
## @ $isMultiBaseNonFrameShiftDel
## @ makeKMerList($s,$k,$i,$j)
sub makeKMerList
{
	my ($s,$k,$i,$j,$bgPepList,$isMultiBaseNonFrameShiftDel)=@_;
	if(!defined($isMultiBaseNonFrameShiftDel)) { $isMultiBaseNonFrameShiftDel=0; }
	my @s=split //,$s;
	my @ret=();
	for(my $t=$i-($k-1);$t<=($j-$isMultiBaseNonFrameShiftDel) && $t<=@s-($k-1);$t++)
	{
		if($t<1) { next; }
		my $m="";
		if($i==$j) { $m=$i-$t+1; }
		else 
		{
		  	my $p1=$i-$t+1;
			my $p2=$j-$t+1;
			if($p2>$k) 
			{ 
				$p2=$k;
			}
			if($p1<1)
			{
				$p1=1;
			}
			if($p1==$p2) { $m="$p1"; }
			else { $m="$p1..$p2"; }
		}
		my $pep=join('',@s[$t-1 .. $t+$k-2]);
		if(!$bgPepList->{$pep})
		{
			push @ret,[$pep,$m];
		}
		
	}
	return @ret;
}

## @ get all possible (2^n) combination of two array. With a test code commented
sub getAllPossibleCombination
{
	my ($pAry1,$pAry2,$n)=@_;
	if($n==1)
	{
		return ($pAry1,$pAry2);
	}elsif($n>1)
	{
		my @_a1=@$pAry1[1 .. $n-1];
		my @_a2=@$pAry2[1 .. $n-1];
		my @_c=getAllPossibleCombination(\@_a1,\@_a2,$n-1);
		my @ret=();
		foreach (@_c)
		{
			push @ret,[$pAry1->[0],@{$_}];
			push @ret,[$pAry2->[0],@{$_}];
		}
		return @ret;
	}
}

sub formatMutType
{
	my ($mutType)=@_;
	if($mutType==SINGLEBASE_SUBSTITUTION) { return "SINGLEBASE_SUBSTITUTION"; }
	elsif($mutType==SINGLEBASE_DELETION) { return "SINGLEBASE_DELETION"; }
	elsif($mutType==MULTIBASE_DELETION) { return "MULTIBASE_DELETION"; }
	elsif($mutType==BLOCK_SUBSTITUTION) { return "BLOCK_SUBSTITUTION"; }
	elsif($mutType==INSERTION) { return "INSERTION"; }
}


sub translateDNA {
	my ($seq) = @_;
	my ($nt3, $protein);
	$seq = uc $seq;
	#length ($seq) % 3 == 0 or printerr "WARNING: length of DNA sequence to be translated is not multiples of 3: <length=${\(length $seq)}>\n";
	while ($seq =~ m/(...)/g) {
		defined $codon1{$1} or print "WARNING: invalid triplets found in DNA sequence to be translated: <$1> in <$seq>\n" and die;
		$protein .= $codon1{$1};
	}
	return $protein;
}

## @ mutate the @$pAryDNA (and also @$pAryDNAW, @$pAryDNAWDis if defined)
## @ $start,$end,$ref,$obs define one mutation
## @ return aastart (aa position of variant start point in original protein;
## @        aaend1  (aa position of variant end point in original protein;
## @        aaend2  (aa position of variant end point in mutant protein;
sub mutateSequenceArray
{
	my ($start,$end,$ref,$obs,$pAryDNA,$pAryDNAW,$pAryDNAWDis)=@_;

	my $aastart = int(($start-1)/3)+1;
	my $aaend1 = int(($end-1)/3)+1;
	my $aaend2;
	my $mutType;

	if ($start == $end and not $ref and $obs) {		#this is an insertion
		$pAryDNA->[$start-1].=$obs;
		$pAryDNAWDis->[$start-1].="-" x length($obs);
		$aaend2=int(($start+length($obs)-1)/3)+1;
		$mutType=INSERTION;
	}elsif($start == $end and $ref and $obs) {		#this is an single base substitution
		$pAryDNA->[$start-1]=$obs;
		$aaend2=int(($start-1)/3)+1;
		$mutType=SINGLEBASE_SUBSTITUTION;
	}elsif($start == $end and $obs eq '') {			#this is an single base deletion
		$pAryDNA->[$start-1]=$obs;
		$aaend2=int(($start-1)/3)+1;
		$mutType=SINGLEBASE_DELETION;
	}elsif($start != $end and $obs eq '') { 		#this is multi-base deletion
		for my $_j ($start .. $end)
		{
			$pAryDNA->[$_j-1]='';
		}
		$aaend2=int(($start-1)/3)+1;
		$mutType=MULTIBASE_DELETION;
		$pAryDNAW->[$start-1]=join('',@$pAryDNAW[$start-1 .. $end-1]);
		$pAryDNAWDis->[$start-1]=join('',@$pAryDNAWDis[$start-1 .. $end-1]);
		for my $_j ($start+1 .. $end)
		{
			$pAryDNAW->[$_j-1]='';
			$pAryDNAWDis->[$_j-1]='';
		}
	}elsif($start != $end and $obs ne '') {			#this is a block substitution
		$pAryDNA->[$start-1]=$obs;
		for my $_j ($start+1 .. $end)
		{
			$pAryDNA->[$_j-1]='';
		}
		$aaend2=int(($start+length($obs)-1-1)/3)+1;
		$mutType=BLOCK_SUBSTITUTION;
		$pAryDNAW->[$start-1]=join('',@$pAryDNAW[$start-1 .. $end-1]);
		$pAryDNAWDis->[$start-1]=join('',@$pAryDNAWDis[$start-1 .. $end-1]);
		for my $_j ($start+1 .. $end)
		{
			$pAryDNAW->[$_j-1]='';
			$pAryDNAWDis->[$_j-1]='';
		}
	}
	return ($aastart,$aaend1,$aaend2,$mutType);
}


=head1 SYNOPSIS

 coding_change.pl [arguments] <exonic-variant-function-file> <gene-def-file> <fasta-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
        -p, --protein=<file>            output all proteins to <file>
        -k, --kmer                      setting kmer [default 9]
        -e, --peptide=<file>            output kmer(peptide) to <file>
        -i, --independent               variants in the same transcript are independent [default NOT]
        -s, --sample                    sample id [default ""]
        --vcf                           <exonic-variant-function-file> is in VCF format [default NOT]
	


 Function: infer the translated protein sequence (or mRNA sequence) for exonic 
 frameshift mutations (or SNPs) identified by ANNOVAR
 
 Example: coding_change.pl ex1.avinput.exonic_variant_function humandb/hg19_refGene.txt humandb/hg19_refGeneMrna.fa
 
 Version: $Date: 2015-03-22 15:29:51 -0700 (Sun, 22 Mar 2015) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--vcf>

<exonic-variant-function-file> is in VCF format.

=item B<--man>

print the complete manual of the program.

=item B<--protein=<file>>

output all proteins to <file>

=item B<--peptide=<file>>

output kmer(peptide) to <file>

=item B<--kmer=<file>>

seting kmer length(k)

=item B<--verbose>

use verbose output.

=item B<--sample>

sample id.

=item B<--mrnaseq>

print out mRNA sequences rather than protein sequences.

=back

=head1 DESCRIPTION

This program will infer the protein sequence for frameshift mutations identified 
by ANNOVAR. Typically, for non-synonymous mutations, the annotate_variation.pl 
program in ANNOVAR will report the amino acid change at the position; however, 
for frameshift mutations which may affect a long stretch of amino acid, ANNOVAR 
will only give a frameshift annotation without printing out the new protein 
sequence. This program can take ANNOVAR exonic_variant_function file and attempt 
to infer the new protein sequence.

=over 8

=item * B<Known Bug>

This program does not handle mitochondria mutations correctly yet.

=back

For questions, comments or bug reports, please contact me at 
$Author: Kai Wang <kai@openbioinformatics.org> $.

=cut
