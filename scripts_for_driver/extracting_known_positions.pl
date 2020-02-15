#!/usr/bin/perl

use strict;
use warnings;

#This script takes the output of calling_majorallele.pl
#and extracts only snp positions present in a known variation vcf file
#Author: Kara (early phd don't laugh, someday I'll come back and clean this up)

chomp (my $IN1 = $ARGV[0]); #known variations (standards)
chomp (my $IN2 = $ARGV[1]); #any vcf file
chomp (my $OUT = $ARGV[2]); #output file (vcf file with only sites from known variations

open (INFILE, "<", $IN1) || die "Cannot open known variations list: $!\n";
my %snps;
for my $line1 (<INFILE>) {
	chomp $line1;
	next if $line1 =~ /^\s*#/;
	my @line1 = split (/\t/, $line1);
	my $chr1 = shift (@line1);
	my $pos1 = shift (@line1);
	my $whole1 = "$chr1".":"."$pos1";
	if ($chr1 eq "Chr") { next;
	} else {$snps{$whole1} = 1;
#		print "$whole1\n";
	}
}
close INFILE;
print "loaded known variants\n";

open (INFILE2, "<", $IN2) || die "Cannot open known variations file: $!\n";
open (OUT, ">", $OUT) || die "creation failed:$!\n";

for my $line2 (<INFILE2>) {
	chomp $line2;
	if ($line2 =~ /^\s*#/) {
	    print OUT "$line2\n";
	} else {

	    my @line2 = split(/\t/, $line2);
	    my $chr2 = shift (@line2);
	    my $pos2 = shift (@line2);
	    my $whole2 = "$chr2".":"."$pos2";
	    if ($chr2 eq "Chr") {next;}  
          if (defined $snps{$whole2}) {
		print OUT "$line2\n";
	    }
	}
}

close INFILE2;
close OUT;

print "Read through entire script\n";


