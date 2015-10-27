#!/usr/bin/perl -w 

=pod

=head1 NAME

			reformat_annovar.pl

=head1 SYNOPSIS

perl reformat_annovar.pl <infile> -id [ID]

=head1 DESCRIPTION

reformat table_annovar output file to pretty table format

=head1 OPTIONS

=over

=item -id

specify the ID of the samples, must be in the same order with annovar input

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($id,$h,$m); 
GetOptions(
	"id=s"=>\$id,
	"h|help"=>\$h,
	"m|man"=>\$m
);

if ($h) {pod2usage(-verbose => 1);}
if ($m) {pod2usage(-verbose  => 2);}
if (!defined $h and !defined $m and !defined $id)  { pod2usage(-verbose=>0);}


my $infile = shift;
my @id = split ",", $id;

open In, "$infile" or die "$!";
my @fields = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER");
my ($length,%header); 

while (<In>) {
	chomp;
	if (/^Chr/) {
		my @new_fields = split "\t", $_;
		$length = @new_fields;
		splice (@new_fields, 0, 5);
		pop @new_fields;
		@fields = (@fields, @new_fields,"INFO","FORMAT",@id);
		my $header = join "\t", @fields;
		print "$header\n";
	}
	else {
		my @records = split "\t", $_;
		print "$records[$length]";
		foreach (1..6) {
			print "\t$records[$length+$_]";
		}
#		print "\t$records[$length-1]";
		if ($records[5] eq "splicing") {
			if ($records[6] =~ /\((.+?)\)/) {
				$records[8] = $1;
				$records[6] =~ s/\(.+?\)//;
			}
		}	
		foreach (5..$length-2) {
			print "\t$records[$_]";
		}
		print "\t$records[$length+7]";
		foreach (($length+8)..$#records) {
			print "\t$records[$_]";
		}
		print "\n";
	}
}
close In;
