#!/usr/bin/perl

# Septeber 29, 2017

#Author: Kara

# ABOUT SCRIPT
#  Script is based off the original callingmajoralle.pl script
#  Input: a vcf file with het calls
#  Output: a vcf file with het calls assined as hom. calls
#
#  Updates:
#  -Fixed genotypes with 1/2 which were not fixed with old version
#  -Now looks at percentage of allele depth and assigns a hom. call
#   if the % of an allele is over 60% (old version just took 
#   whichever one was higher

use strict;
use warnings;

chomp (my $IN1 = $ARGV[0]); #any vcf file with het. calls
chomp (my $OUTFILE = $ARGV[1]); #output vcf files with corrected het calls (major allele called)

open (INFILE, "<", $IN1) || die "Cannot open known variations list: $!\n";
open (OUTFILE, ">", $OUTFILE) || die "Cannot open output file:$!\n";

#Reading in vcf file with het calls

for my $line1 (<INFILE>) {
	chomp $line1;
	
#Assigning variables (def a more efficient way to do this)

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

#If a genotype is a het call (0/1 or 1/2), then look at the frequencies
#Assign a new homozygous genotype call with the allele that has a read depth >= 60%
#If a site does not have one allele with an allele >= 60%,
#  then assign a missing (./.) genotype

	    my $new_geno;	
	    if ($geno eq "0/1") {
		my ($freq1, $freq2) = split /,/, $ad;
		my $perc1 = $freq1 / ($freq1 + $freq2);
		my $perc2 = $freq2 / ($freq1 + $freq2);
		if ($perc1 >= 0.70) {
		    $new_geno = "0/0";
		} elsif ($perc2 >= 0.70) {
		    $new_geno = "1/1";
		} else {
		    $new_geno = "./.";
		}
		my $new_sample = $new_geno . ":" . $ad . ":" . $dp . ":" .  $s1 . ":" . $s2;
		print OUTFILE "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\t$info\t$format\t$new_sample\n";

	    } elsif ($geno eq "1/2") {
		my ($freq1, $freq2) = split /,/, $ad;
                my $perc1 = $freq1 / ($freq1 + $freq2);
                my $perc2 = $freq2 / ($freq1 + $freq2);
		if ($perc1 >= 0.70) {
                    $new_geno = "1/1";
		} elsif ($perc2 >= 0.70) {	
		    $new_geno = "2/2";
		} else {
		    $new_geno = "./.";
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
