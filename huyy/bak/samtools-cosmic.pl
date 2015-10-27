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
 
open IN, "/PROJ/GR/share/medinfo.00database/annovar/humandb_b37/hg19_cosmic65.txt" or die "ERROR: open $pesoap: $!";
if($vcf=~/gz$/) {
	open INN,"/PROJ/GR/share/medinfo.02pipeline/bin/bgzip -dc $vcf |" or die "ERROR: open $vcf: $!";
}else{
  open INN, $vcf or die "ERROR: open $vcf: $!";
}
open OUTM, ">$outfile" or die "ERROR: $!";
while(defined($line=<IN>)) {
	chomp $line;
	@colu=split(/\t/,$line);
	#1       982818  982819  GG      AA      ID=COSM143736;OCCURENCE=1(skin)
	#1       979279  979279  C       -       ID=COSM1345050;OCCURENCE=1(large_intestine)
	if($colu[3] eq "-"){
		$pos=join ":",$colu[0],$colu[1];
		$hash{$pos}=$colu[5];
		$phen{$pos}=$colu[4];
	}elsif($colu[4] eq "-"){
		$pos=join ":",$colu[0],$colu[1];
		$hash{$pos}=$colu[5];
		$phen{$pos}=$colu[3];
	}
}
close IN;
while(defined($line=<INN>)){
	chomp $line;
	if($line!~/^#/) {
		$mark=0;
		@colu=split(/\t/,$line);
		$st=$colu[1];
		$chr=$colu[0];
				if(length($colu[4])>length($colu[3])) {
			#$type="insertion";
			$sub=$colu[4];
			$quer=$colu[3];
		}else{
			#$type="deletion";
			$sub=$colu[3];
			$quer=$colu[4];
		}
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
				@quer=split(//,$quer);
			foreach (-1,0..$#quer) {
				if($_ == -1) {
					$seq=join("",$phen{$pos},@quer);
				}elsif($ _ == $#quer) {
					$seq=join("",@quer,$phen{$pos});
				}else{
					$seq=join("",@quer[0..$_],$phen{$pos},@quer[$_+1..$#quer]);
				}
				if($seq eq $sub) {
				$mark=1;
				$rate=$hash{$pos};
				last;
			}
			}
		}
	}
		if($mark==0) {
			print OUTM $line,"\n";
		}else{
			if($colu[7]!~/cosmic/){
				$info=$colu[7].";cosmic=".$rate;
			  print OUTM join ("\t",@colu[0..6],$info,@colu[8..$#colu]),"\n";
			}else{
				print OUTM $line,"\n";
			}
			
		}
	}else{
		print OUTM $line,"\n";
	}
}
		
