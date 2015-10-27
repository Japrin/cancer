#!/usr/bin/perl
my $usage=<<INFO;
Descriptions:

	To check if the samples have been contaminated or not
	This perl will generate run_blast_nt.sh at the output dir
	
Command-line Option:
	-f 			the input list of fq, fq.gz, fa, or fa.gz files, each line
				has three columns: fq1\tfq2\tkeyName, the fqFiles should be raw data
	-o			the output dir for the results, default: ./
	-bp			the path of blast, default:/ngs/self-software/blastall	
	-bd	 		the database of blast, default:/ngs/pub/genome/db/blast/nt.
	-nb			the number of database, default:10, means nt.00~nt.nb
	-a			number of processors to use, default:6
	-n			number of reads to run with nt database(n/2 for each fqFile in paired end)
	-d			minimal identity to filter blast results, default:0.95
	-l			minimal mapped length to filter blast results, default:99
	-p			minimal percentage for output results, default:0.005		
	-h			output help information			
	
Usage:

	perl shell4nt.pl -h
	perl shell4nt.pl -f in.list -o blast -i 0.95 -l 99

Version:
	2011.11.21

Contact:

 	Zhang Jinbo Email: zhangjinbo\@novogene.cn
 	
INFO

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use Cwd;
use PerlIO::gzip;
use Cwd 'abs_path';

my ($listFile, $outDir, $blastPath, $blastDatabase, $numberOfDatabase);
my ($processNumber, $numberOfReads, $minIdentity, $minMappedLength, $minPercentage);
my ($help);

GetOptions(
	"f=s"=>\$listFile,
	"o=s"=>\$outDir,
	"bp=s"=>\$blastPath,
	"bd=s"=>\$blastDatabase,
	"a=s"=>\$processNumber,
	"n=i"=>\$numberOfReads,
	"d=f"=>\$minIdentity,
	"l=s"=>\$minMappedLength,
	"p=f"=>\$minPercentage,
	"nb=i"=>\$numberOfDatabase,
	"h"=>\$help,
);

die "$usage" if($listFile eq "");
die "$usage" if(defined $help);

$listFile = abs_path($listFile);
$outDir ||= "./";
$outDir = abs_path($outDir);
$blastPath ||= "/ngs/self-software/blastall";
$blastDatabase ||= "/ngs/pub/genome/db/blast/nt.";
$processNumber ||= 6;
$minIdentity ||= 0.95;
$minMappedLength ||= 99;
$minPercentage ||= 0.005;
$numberOfDatabase ||= 10;
$numberOfReads ||= 10000;

# dir for programe
my $blast_parse = "perl /ngs/PG/zhangjinbo/bin/common/blast_nt/blast_parser.pl";
my $blast_stats = "perl /ngs/PG/zhangjinbo/bin/common/blast_nt/blast_stats.pl";
my $blast_merge = "perl /ngs/PG/zhangjinbo/bin/common/blast_nt/merge4blast_stats.pl";

my %samples = ();

&readFile($listFile, \%samples);


sub readFile{
	my ($inFile, $sub) = @_;
	my $shell = $outDir."/run_blast_nt.sh";
	open OUT, ">".$shell or die $!;
	my %merge = ();
	my $sample_id = 0;
        print $inFile."\n";
	open IN, $inFile or die $!;
	while(<IN>){
		chomp;
		
		$sample_id++;
		
		my @arr = split('\s+', $_);
		my $key = $arr[2];
		
		print STDERR "$listFile:$arr[0] does not exist!\n" if(! -e $arr[0]);
		print STDERR "$listFile:$arr[1] does not exist!\n" if(! -e $arr[1]);
		my $sampleDir = $outDir."/".$key;
		system("mkdir $sampleDir\n") if(! -d $sampleDir);
		# input file
		my $fqDir = $sampleDir."/fq";
		system("mkdir $fqDir\n") if(! -d $fqDir);
		my $in = $fqDir."/".$key.".fq";
		# blast output
		my $blastDir = $sampleDir."/blast";
		system("mkdir $blastDir\n") if(! -d $blastDir);
		# result output
		my $resultDir = $sampleDir."/result";
		system("mkdir $resultDir\n") if(! -d $resultDir);
		
		
		&geneFq($in, $arr[0], $arr[1]);
		# run for blast
		for(my $i = 0; $i <= $numberOfDatabase; $i++){
			my $index = "";
			if($i < 10){
				$index = "0".$i;
			}else{
				$index = $i;
			}
			my $tmp4blastOut1 = $blastDir."/".$key.".".$index;
			my $tmp4blastOut2 = $tmp4blastOut1.".table";
			my $tmp4result = $resultDir."/".$key.".".$index."_db.xls";
			print OUT $blastPath." -d $blastDatabase".$index." -i $in -o $tmp4blastOut1 -a $processNumber -p blastn -e 1e-05  -F F\n";
			print OUT $blast_parse." $tmp4blastOut1 > $tmp4blastOut2\n";
			print OUT $blast_stats." -i $tmp4blastOut2 -o $tmp4result -nb $numberOfReads -mp $minPercentage -mi $minIdentity -ml $minMappedLength\n";
			
		}
		
		my $mergeSH = $resultDir."/merge.sh";
		#print $mergeSH."\n";
		open TMP, ">".$mergeSH or die $!;
		print TMP $blast_merge." ".$resultDir." ".$key."\n";
		close(TMP);
		$merge{"sh ".$mergeSH} = $sample_id;
	}
	close(IN);
	
	my $qsub = $outDir."/qsub_blast.sh";
	open QSUB, ">".$qsub or die $!;
	print QSUB "perl /ngs/PG/zhangjinbo/bin/common/qsub-sge.pl --resource vf=2G -lines 3 $shell --convert no\n";
	# for the merge
	foreach my $key ( sort {$merge{$a} <=> $merge{$b}} keys %merge ){
		print QSUB $key."\n";
	}
	close(QSUB);
	
	print "Just type the following to finish the task:\n";
	print "nohup sh $qsub > $qsub.log 2> $qsub.error &\n";
	
	
}

sub geneFq{
	my($out, $in1, $in2) = @_;
	my $tmp4num = int($numberOfReads/2);
	open FQ, ">".$out or die $!;
	if($in1 =~ /gz$/){
		open IN1, "<:gzip", $in1 or die $!;
		open IN2, "<:gzip", $in2 or die $!;
	}else{
		open IN1, $in1 or die $!;
		open IN2, $in2 or die $!;
	}
	my $tmp4recod = 0;
	while(<IN1>){
		$tmp4recod++;
		my $id1 = $_;
		my $id2 = <IN2>;
		my $seq1 = <IN1>;
		my $seq2 = <IN2>;
		
		print FQ ">".$tmp4recod."_1\n";
		print FQ $seq1;
		print FQ ">".$tmp4recod."_2\n";
		print FQ $seq2;
		
		if($id1 =~ /^>/ && $id2 =~ /^>/){
			# for fa file
			#next;
			
		}elsif($id1 =~ /^@/ && $id2 =~ /^@/){
			# for fq file
			<IN1>;
			<IN1>;
			<IN2>;
			<IN2>;
			
		}else{
			print STDERR "Errors happen at input file : $in1 and $in2\n";
		}
		last if($tmp4recod >= $tmp4num);
	}
	close(IN1);
	close(IN2);
	close(FQ);
}




















