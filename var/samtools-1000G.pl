#!/usr/bin/perl -w
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=head1 Function
  Indel annotation with 1000G
  
=head1 Usage
  perl -i <vcf> -o <outile>
  
=head1 Options
  -i<fasta>         input file in vcf format
  -o<outfile>       output file 
=head1 Author
  Kong gunayi; kongguanyi@novogene.cn
  
=head1 Version
  v1.0; 2013-9-6
  
=cut

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
use Getopt::Long;

my ($pesoap,$out_file);
GetOptions(
  "i=s" =>\$vcf, 
  "o|out-file=s"    =>\$outfile
  );
 die `pod2text $0` if (!$vcf || !$outfile);
 
open IN, "/DBS/DB_temp/zhangLab/annovar/humandb_b37/hg19_ALL.sites.2012_04.txt" or die "ERROR: open $pesoap: $!";
if($vcf=~/gz$/) {
	open INN,"bgzip -dc $vcf |" or die "ERROR: open $vcf: $!";
}else{
  open INN, $vcf or die "ERROR: open $vcf: $!";
}
open OUTM, ">$outfile" or die "ERROR: $!";
while(defined($line=<IN>)) {
	chomp $line;
	@colu=split(/\t/,$line);
	$colu[0]=~s/^chr//;
	if($colu[3]=~/^\d/){
	$pos=join ":",$colu[0],$colu[1];
	$hash{$pos}=$colu[4];
}
}
close IN;
while(defined($line=<INN>)){
	chomp $line;
	if($line!~/^#/) {
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
		#if($opt_a)
		{
			$_l=$st;
			$_r=$st+$ll-1;
		}
		#else
		#{
		#	$_l=$pos-$opt_d;
		#	$_r=$pos+$opt_d;
		#}
		#for($i=$st-3;$i<=$st+3;$i++){
		for($i=$_l;$i<=$_r;$i++){
			$pos=join ":",$chr,$i;
			if(defined($hash{$pos})) {
				$mark=1;
				$rate=$hash{$pos};
				last;
			}
		}
		if($mark==0) {
			print OUTM join("\t",@colu),"\n";
		}else{
			if($colu[7]!~/1000g2012apr_all/i){
				$info=$colu[7].";1000g2012apr_all=".$rate;
			  print OUTM join ("\t",@colu[0..6],$info,@colu[8..$#colu]),"\n";
			}else{
				print OUTM join("\t",@colu),"\n";
			}
			
		}
	}else{
		print OUTM $line,"\n";
	}
}
		
