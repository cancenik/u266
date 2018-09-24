#!/usr/bin/perl -w 

## CC 
################
# This script generates the CDS count table
## We will also output coverage and length of the region
################

use strict; 
use File::Basename;
use constant ID => 4; 
use constant COUNT => 6; 
use constant COVERED_BASES => 7; 
use constant FEATURE_LEN => 8; 
use constant REGION => 3; 
use Data::Dumper; 
$Data::Dumper::Terse = 1;
$Data::Dumper::Indent = 0;
require Simple_Utils;

my %counts;
my @files = </srv/gsfs0/projects/snyder/ccenik/RIBOSEQ_OTHERS/U266_Myeloma/*/*/Transcript_Counts*>; 

my $counts_file = "Transcript_Counts_All_Libraries.tsv"; 
unless (-e $counts_file) {
open (OUT, ">$counts_file") || die "Cannot open output file: $!\n"; 
print OUT "ID\tREGION\t"; 
foreach my $file (@files) { 
    my $dirname = dirname($file);
    my $basename = basename($dirname);
    print OUT "$basename"."_Counts\t$basename" ."_CoveredBases\t$basename". "_FeatureLen\t"; 
#    open (IN, "gzip -dc $file |"); 
    open (IN, "$file") || die "Cannot open file: $!\n"; 
    while (<IN>) { 
	my @F = split; 
	my @ID = split (/\|/, $F[0]); 
	push (@{$counts{$ID[ID]}->{$F[REGION]}}, ($F[COUNT], $F[COVERED_BASES], $F[FEATURE_LEN])); 
    }
}
print OUT "\n"; 
#print OUT Dumper(\%counts); 

my $out_fh = *OUT;
Simple_Utils::print_hash(\%counts, 0 , $out_fh);
close OUT; 

#REFORMAT OUTPUT FOR R_INPUT
open (REFORMAT, $counts_file) || die "Cannot open file for reformatting\n"; 
my $header = <REFORMAT>; 
my $reformatted_file = "Reformatted_Transcript_Counts_All_Libraries.tsv";
#my $reformatted_file = "Reformatted_Species_Counts_All_Libraries.tsv";
open (FINALOUT, ">$reformatted_file") || die "Cannot open final output: $!\n";
print FINALOUT $header;
my $tmp_id; 
while (<REFORMAT>) {
    chomp;
    my @G = split (/\t/); 
    if ($G[0] eq "") {
	$G[2] =~ s/\s+/\t/g;
	print FINALOUT "$tmp_id\t$G[1]\t$G[2]\n";
    }
    else {
	$tmp_id = $G[0]; 
    }
}
close REFORMAT;
close FINALOUT; 
}


