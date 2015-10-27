#!/usr/bin/perl -w
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head1 Function
  Indel annotation with dbsnp
  
=head1 Usage
  perl -i <vcf> -o <outile>
  
=head1 Options
  -i<fasta>         input file in vcf format
  -o<outfile>       output file 
  -m			mode, SNP or INDEL [default INDEL]
=head1 Author
  Kong gunayi; kongguanyi@novogene.cn
  
=head1 Version
  v1.0; 2013-9-6
  
=cut

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use Getopt::Long;

my ($pesoap,$out_file,$opt_m);
GetOptions(
  "i=s" =>\$vcf, 
  "o|out-file=s"    =>\$outfile,
  "m=s" =>\$opt_m
  );
 die `pod2text $0` if (!$vcf || !$outfile);
if(!defined($opt_m)) { $opt_m="INDEL"; } 
 
open IN, "/DBS/DB_temp/zhangLab/annovar/humandb_b37/hg19_snp138.txt" or die "ERROR: open $pesoap: $!";
if($vcf=~/gz$/) {
	open INN,"bgzip -dc $vcf |" or die "ERROR: open $vcf: $!";
}else{
  open INN, $vcf or die "ERROR: open $vcf: $!";
}
open OUTM, ">$outfile" or die "ERROR: $!";
while(defined($line=<IN>)) {
	chomp $line;
	@colu=split(/\t/,$line);
	$colu[1]=~s/^chr//;
	if($opt_m eq "INDEL")
	{
		if($colu[11]!~/single/)
		{
			$pos=join ":",$colu[1],$colu[2];
			$hash{$pos}=$colu[11];
			$rs{$pos}=$colu[4];
		}
	}elsif($opt_m eq "SNP")
	{
		if($colu[11]=~/single/)
		{
			$pos=join ":",$colu[1],$colu[3];
			$hash{$pos}=$colu[11];
			$rs{$pos}=$colu[4];
		}
	}
}
close IN;
while(defined($line=<INN>))
{
	chomp $line;
	if($line!~/^#/) 
	{
		$mark=0;
		@colu=split(/\t/,$line);
		$st=$colu[1];
		$colu[0]=~s/^chr//;
		$chr=$colu[0];
		
		my $ref=$colu[3];
		my $alt=$colu[4];
		my $len_m=length($ref);
		my $len_n=length($alt);
		my $ll=$len_m>$len_n?$len_m:$len_n;
		my ($_l,$_r);
		{
			$_l=$st;
			$_r=$st+$ll-1;
		}
		
		for($i=$_l;$i<=$_r;$i++){
			$pos=join ":",$chr,$i;
			if(defined($hash{$pos})) {
				$mark=1;
				$rs=$rs{$pos};
				last;
			}
		}
		if($mark==0) {
			print OUTM join("\t",@colu),"\n";
		}else{
			print OUTM join ("\t",@colu[0..1],$rs,@colu[3..$#colu]),"\n";
		}
	}else{
		print OUTM $line,"\n";
	}
}
		
