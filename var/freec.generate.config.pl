#!/usr/bin/env perl
#============================================================================
# Name        		: freec.generate.config.pl
# Author      		: zhengliangtao
# Version     		: v1.00
# Created On  		: Tue Dec 24 18:11:00 2013
# Last Modified By	: 
# Last Modified On	: Tue Dec 24 18:11:00 2013
# Copyright   		: Copyright (C) 2013
# Description 		: 
#============================================================================

=pod

=head1 Usage

	perl freec.generate.config.pl [option] <tumor bam> [<normal bam>]
	
	-o			outDir [default . ]
	-m			mateOrientation, 0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends) [default FR]
	-e			example config file [default "/PROJ/GR/share/medinfo.02pipeline/cancer/var/freec.config.example" ]
	-gemMappabilityFile	mappability file [default "/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/gem/human_g1k_v37_decoy.mappability.100bp.out.mappability" ]
	-chrFiles		chrFiles [default "/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/byChr/" ]
	-chrLenFile		chrLenFile [default "/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/chr.24.length" ]
	-gender			gender [default "XX"]
	-contaminationAdjustment	contaminationAdjustment [default "TRUE" ]
	-SNPfile		SNPfile [default "/PROJ/GR/share/medinfo.00database/testData/freec/b37_snp137.SingleDiNucl.1based.txt" ]
	-h			display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Path;

my ($in,$out);
my ($opt_h,$opt_o,$opt_e,$opt_gem,$opt_chrFiles,$opt_chrLenFile,$opt_gender,$opt_contaminationAdjustment,$opt_SNPfile,$opt_m);
GetOptions("h"	=>\$opt_h,"o=s"=>\$opt_o,"e=s"=>\$opt_e,"gemMappabilityFile=s"=>\$opt_gem,"chrFiles=s"=>\$opt_chrFiles,"chrLenFile=s"=>\$opt_chrLenFile,"gender=s"=>\$opt_gender,"contaminationAdjustment=s"=>\$opt_contaminationAdjustment,"SNPfile=s"=>\$opt_SNPfile,"m=s"=>\$opt_m);
if(@ARGV<1 || $opt_h) { usage(); }
if(!defined($opt_o)) { $opt_o="."; }
if(!defined($opt_e)) { $opt_e="/PROJ/GR/share/medinfo.02pipeline/cancer/var/freec.config.example"; }
if(!defined($opt_gem)) { $opt_gem="/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/gem/human_g1k_v37_decoy.mappability.100bp.out.mappability"; }
if(!defined($opt_chrFiles)) { $opt_chrFiles="/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/byChr/"; }
if(!defined($opt_chrLenFile)) { $opt_chrLenFile="/PROJ/GR/share/medinfo.00database/genome/human/b37_gatk/chr.24.length"; }
if(!defined($opt_gender)) { $opt_gender="XX"; }
if(!defined($opt_contaminationAdjustment)) { $opt_contaminationAdjustment="TRUE"; }
if(!defined($opt_SNPfile)) { $opt_SNPfile="/PROJ/GR/share/medinfo.00database/testData/freec/b37_snp137.SingleDiNucl.1based.txt"; }
if(!defined($opt_m)) { $opt_m="FR"; }

mkpath $opt_o;

open $in,$opt_e or die "Cann't open file $opt_e ($!) \n";
#open $out,">","$opt_o/freec.config" or die "$!";
while(<$in>)
{
	chomp;
	my $line=$_;
	if(/^\s*$/ || /^#/) { print "$line\n"; next; }
	if(/^outputDir/) { print "outputDir = $opt_o\n"; }
	elsif(/^gemMappabilityFile/) { print "gemMappabilityFile = $opt_gem\n"; }
	elsif(/^chrFiles/) { print "chrFiles = $opt_chrFiles\n"; }
	elsif(/^chrLenFile/) { print "chrLenFile = $opt_chrLenFile\n"; }
	elsif(/^sex/) { print "sex = $opt_gender\n"; }
	elsif(/^contaminationAdjustment/) { print "contaminationAdjustment = $opt_contaminationAdjustment\n"; }
	elsif(/^SNPfile/) { print "SNPfile = $opt_SNPfile\n"; }
	elsif(/^\[sample\]/)
	{
		print "[sample]\n";
		print "mateFile = $ARGV[0]\n";
		##print "inputFormat = BAM\n";
		print "mateOrientation = $opt_m\n";
	}elsif(/^\[control\]/)
	{
		if($ARGV[1])
		{
			print "[control]\n";
			print "mateFile = $ARGV[1]\n";
			##print "inputFormat = BAM\n";
			print "mateOrientation = $opt_m\n";
		}
	}else
	{
		print "$line\n";
	}

}
###################################################################




sub usage
{
	die `pod2text $0`;
}
