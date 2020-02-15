#!/usr/bin/perl

use strict;
use warnings;

chomp (my $IN1 = $ARGV[0]); #any vcf file with het. calls
chomp (my $OUTFILE = $ARGV[1]); #output vcf files with corrected het calls (major allele called)

open (INFILE, "<", $IN1) || die "Cannot open known variations list: $!\n";
open (OUTFILE, ">", $OUTFILE) || die "Cannot open output file:$!\n";

for my $line1 (<INFILE>) {
	chomp $line1;
	
	if ($line1 =~ /^\s*#/) {
	    print OUTFILE "$line1\n";
	} else {
	    my @line1 = split (/\t/, $line1);
	    my $chrom = shift @line1;
	    my $pos = shift @line1;
	    my $id = shift @line1;
	    my $ref = shift @line1;
	    my $alt = shift @line1;
	    my $qual = shift @line1;
	    my $filter = shift @line1;
	    my $info = shift @line1;
	    my $format = shift @line1;
	    my $sample = shift @line1;
	    my ($geno, $ad, $dp, $s1, $s2) = split /:/, $sample;
	    #print "$ref\t$alt\t$info\n";

	    my $new_geno;	
	    if ($geno eq "0/1") {
		my ($freq1, $freq2) = split /,/, $ad;
		if ($freq1 >= $freq2) {
		    $new_geno = "0/0";
		} else {
		    $new_geno = "1/1";
		}
		my $new_sample = $new_geno . ":" . $ad . ":" . $dp . ":" .  $s1 . ":" . $s2;
		print OUTFILE "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$new_sample\n";
	    } else {
		print OUTFILE "$line1\n";
	    }
	}
}

close $IN1;
close $OUTFILE;

print "Read through entire script\n";
