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
 
open IN, "/PROJ/GR/share/medinfo.00database/annovar/humandb_b37/hg19_snp137.txt" or die "ERROR: open $pesoap: $!";
if($vcf=~/gz$/) {
	open INN,"/PROJ/GR/share/medinfo.02pipeline/bin/bgzip -dc $vcf |" or die "ERROR: open $vcf: $!";
}else{
  open INN, $vcf or die "ERROR: open $vcf: $!";
}
open OUTM, ">$outfile" or die "ERROR: $!";
while(defined($line=<IN>)) {
	chomp $line;
	@colu=split(/\t/,$line);
	#if($colu[11]=~/deletion|insertion/){
	if($opt_m eq "INDEL")
	{
		if($colu[11]!~/single/)
		{
			$pos=join ":",$colu[1],$colu[2];
			$hash{$pos}=$colu[11];
			$rs{$pos}=$colu[4];
			@phen=split(/\//,$colu[9]);
			$phen{$pos}=$phen[1];
		}
	}elsif($opt_m eq "SNP")
	{
		if($colu[11]=~/single/)
		{
			$pos=join ":",$colu[1],$colu[3];
			$hash{$pos}=$colu[11];
			$rs{$pos}=$colu[4];
			@phen=split(/\//,$colu[9]);
			$phen{$pos}=$phen[1];
		}
	}
}
close IN;
while(defined($line=<INN>)){
	chomp $line;
	if($line!~/^#/ && length($colu[4])!= length($colu[3])) {
		$mark=0;
		@colu=split(/\t/,$line);
		$st=$colu[1];
		$chr="chr".$colu[0];
		if(length($colu[4])>length($colu[3])) {
			$type="insertion";
			$sub=$colu[4];
			$quer=$colu[3];
		}else{
			$type="deletion";
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
			if(defined($hash{$pos}) && $hash{$pos} eq $type) {
			#if(defined($hash{$pos})) {
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
				$rs=$rs{$pos};
				last;
			}
		}
	}
		}
		if($mark==0) {
			print OUTM $line,"\n";
		}else{
			print OUTM join ("\t",@colu[0..1],$rs,@colu[3..$#colu]),"\n";
		}
	}else{
		print OUTM $line,"\n";
	}
}
		
