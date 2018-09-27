#!/usr/bin/env perl
#============================================================================
# Name        		: addAnn.pl
# Author      		: zhengliangtao
# Version     		: v1.1
# Created On  		: Mon Sep 24 16:09:56 2012
# Last Modified By	: huyiyao
# Last Modified On	: Wed Jul 24 11:02:56 2013
# Copyright   		: Copyright (C) 2012
# Description 		: add an extra annotation according to gene name
#============================================================================

=pod

=head1 NAME

addAnn_hyy.pl  - add annotaion information according to a column

=head1 SYNOPSIS 

perl addAnn.pl -annFile [annotation file] -annName [annotation name] -field [annotation according to field] -help -man <infile>

=head1 AUTHOR

Huyiyao
2013-10-23

=head1 OPTIONS

=over 8

=item	-annFile

the location of the annotation databse, [default: /WPS/BP/zhenglt/00.database/DAVIDKnowledgeBase/other/OMIM.txt] 

=item -annName

Specify the name for the annotation [default: OMIM]
other name include:
OMIM_DAVID
CancerGene
BIND
BBID
COG_NAME
COG_ONTOLOGY
GENE_NAME
GO_BP
GO_CC
GO_MF
KEGG_PATHWAY
PANTHER_PATHWAY
REACTOME_PATHWAY
PUBMED_ID
SMART
UCSC_TFBS

=item -field

Specify which field to be the key for sorting. [default: Gene]

=item -help

Show help information

=item -man

show detailed help information 

=back

=cut


use strict;
use warnings;
#use lib qw(/PROJ/GR/share/Software/medinfo/04lib/perl/lib/site_perl/5.18.0 );
use Getopt::Long;
use Pod::Usage;
use POSIX qw(strftime);


###my $DAVID="/WPS1/zhenglt/00.database/DAVIDKnowledgeBase";
my $DAVID="/WPSnew/zhenglt/00.database/DAVID/DAVIDKnowledgebase";


	my ($in,$out);
	my ($opt_h,$opt_m,$field,$annFile,$annName);
	GetOptions("help"=>\$opt_h,"man"=>\$opt_m,"field=s"=>\$field,"annFile=s"=>\$annFile,"annName=s"=>\$annName);
	
	if(@ARGV<0 || $opt_h) { pod2usage(-verbose=> 1); }
	if ($opt_m) {pod2usage(-verbose=> 2);}
	if (!defined $opt_h and !defined $opt_m and !defined $field and !defined $annFile and !defined $annName ) {pod2usage(-verbose=> 0)};
	if (!defined($annName)) {$annName = "OMIM";}
	if(!defined($field)) { $field="Gene"; }
    ###if (($annName eq "CancerGene") and (!defined $annFile)) {$annFile = "$DAVID/other/CancerGene.txt";}
	if (($annName eq "CancerGeneCensus") and (!defined $annFile)) {$annFile = "$DAVID/";}
	if (($annName eq "OMIM") and (!defined $annFile)) { $annFile = "$DAVID/other/OMIM.txt";}
	if (($annName eq "BIND") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2BIND.txt";}
	if (($annName eq "COG_NAME") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2COG_NAME.txt";}
	if (($annName eq "GENE_NAME") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2DAVID_GENE_NAME.txt";}
	if (($annName eq "PUBMED_ID") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2PUBMED_ID.txt";}
	if (($annName eq "SMART") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2SMART.txt";}
	if (($annName eq "UCSC_TFBS") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2UCSC_TFBS.txt";}
	if (($annName eq "BBID") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2BBID.txt";}
	if (($annName eq "COG_ONTOLOGY") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2COG_ONTOLOGY.txt";}
	if (($annName eq "GO_BP") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2GOTERM_BP_FAT.txt";}
	if (($annName eq "GO_CC") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2GOTERM_CC_FAT.txt";}
	if (($annName eq "GO_MF") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2GOTERM_MF_FAT.txt";}
	if (($annName eq "KEGG_PATHWAY") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2KEGG_PATHWAY.txt";}
	if (($annName eq "OMIM_DAVID") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2OMIM_DISEASE.txt";}
	if (($annName eq "PANTHER_PATHWAY") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2PANTHER_PATHWAY.txt";}
	if (($annName eq "REACTOME_PATHWAY") and (!defined $annFile)) { $annFile = "$DAVID/OFFICIAL_GENE_SYMBOL2REACTOME_PATHWAY.txt";}
	

	my %annList=();
	readAnnFile(\%annList,$annFile);

	my $hLine=<>;
	chomp $hLine;
	print "$hLine\t$annName\n";
	my @Lines = split "\t", $hLine;
	my $col;
	foreach my $i (0..$#Lines) {
		if ($Lines[$i] eq $field) {
			$col = $i;
		}
	}
	if (!defined ($col)) {print "$field not found!\n";}
	
	while(<>)
	{
		chomp;
		my $_line=$_;
		if(/^\s*$/) { next; }
		next if (/#/);
		my @field=split /\t/;
		my $id=$field[$col];
		$id=~s/\(.+?\)//g;
		$id = lc($id);
		my @g=split /[,;]/,$id;
		#printf "%d\t$id\t%s\n",scalar @g,join("\t",@g);
		#next;
		my @ann=();
		my $flag=0;
		foreach (@g)
		{
			if(exists($annList{$_}))
			{
				#push @ann,sprintf "${_}:%s",join(";",@{$annList{$_}});
				push @ann,sprintf "%s",join(";",@{$annList{$_}});
				$flag=1;
			}
			#else
			#{
			#	push @ann,".";
			#}
		}
		printf "%s\t%s\n",$_line,$flag?join(";",@ann):".";

	}


sub readAnnFile
{
	my ($in);
	my ($pList,$infile)=@_;
	open $in,$infile or die "Cann't open file $infile ($!) \n";
	while(<$in>)
	{
		chomp;
		if(/^\s*$/ || /^#/) { next; }
		my @field=split /\t/;
		$field[0] = lc($field[0]);
		if(!exists($pList->{$field[0]})) { $pList->{$field[0]}=[]; }
		push @{$pList->{$field[0]}},$field[1];
	}

}
sub usage
{
	die `pod2text $0`;
}
