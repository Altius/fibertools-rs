#!/bin/env perl
use strict;
use warnings;
use List::Util qw( sum min max );

if (@ARGV != 5) {
    die "Usage:  $0 genome_reference_fasta chromosome_name region_min0 region_max1 output_folder\n";
}
my ($referenceFastaFile, $chromosome, $highlightMin0, $highlightMax1, $outputFolder) = @ARGV;

#$ samtools faidx /home/ehaugen/refseq/hg38/hg38.fa chr11:1010101-1015000 | fastastat.pl
#>chr11:1010101-1015000  4900    bp      4900    non-N bases.

unless (-d $outputFolder) {
    mkdir $outputFolder or die "Failed to create output folder $outputFolder\n";;
}

# Pull the reference just for the region of interest
my $highlightMin1 = $highlightMin0 + 1;
my $coords_one_based = "$chromosome:$highlightMin1-$highlightMax1";
my $referenceString = "";
open (my $infa, "samtools faidx $referenceFastaFile $coords_one_based |") or die "$!";
while (<$infa>) {
    chomp;
    if (/^>(.*)$/) {
        warn "Getting sequence $1\n";
    } else {
        $referenceString .= $_;
    }
}
close $infa;
my $expectedReferenceLength = $highlightMax1 - $highlightMin0;
my $foundReferenceLength = length( $referenceString );
if ($foundReferenceLength != $expectedReferenceLength) {
    die "Expecting $expectedReferenceLength bp but read $foundReferenceLength\n";
}

# Mark the A's and T's that could get m6A calls
my @refBasesAT = ();
for (my $i = 0; $i < $expectedReferenceLength; ++$i) {
    my $base = uc( substr($referenceString, $i, 1 ) );
    if ("A" eq $base or "T" eq $base) {
        push @refBasesAT, $highlightMin0 + $i;
    }
}
my $countATsInRange = scalar(@refBasesAT);

my @hashes = ();
my $lines = 0;
my %vectorOffsets = ();

while (<STDIN>) {
    chomp;
    my $trackLine = $_;
    my ($chrom, $chromStart, $chromEnd, $name, $coverage, $strand, $thickStart, $thickEnd, $itemRgb,
        $blockCount, $blockSizes, $blockStarts ) = split /\t/;
    ++$lines;
    my $mappedLength = $chromEnd - $chromStart;
    # Will include with NA's for missing extent: if ($chromStart > $highlightMin0 or $chromEnd < $highlightMax1) {

    # Not relative locations, just literally keep track of every base
    # and remember to ignore the end two positions of the fiber record
    my @starts = split /\,/, $blockStarts;
    my $ignoreFirst = shift @starts;
    my $ignoreLast = pop @starts;
    next unless (@starts > 0); # no methylation calls at all?  not a useful fiber!
    my %hash_m6A = map { $_ => 1 } @starts;

    my @sortvector = ();
    for my $refBase0 (@refBasesAT) {
        # Ignoring the BED12 end positions that aren't considered in the track
        if ($refBase0 > $chromStart + 1 and $refBase0 < $chromEnd - 1) {
            my $offset = $refBase0 - $chromStart;
            my $methylated = defined( $hash_m6A{$offset} ) ? 1 : 0;
            push @sortvector, $methylated;
        } else {
            # out of range for this molecule, incomplete overlap
            push @sortvector, "NA";
        }
    }

    # Skip this if no A/T bases in range
    next unless (scalar(@sortvector) > 0);

    my $mean = sum(@sortvector) / scalar(@sortvector);
    my %hash = (
        name => $name,
        line => $trackLine,
        vector => \@sortvector,
        inputorder => $lines,
        meanmeth => $mean
    );
    unless (defined($countATsInRange)) {
         $countATsInRange = scalar(@sortvector);
    }
    push @hashes, \%hash;
}

# Matrix of the m6A statuses within excerpt region
my $fileMatrix = "$outputFolder/matrix_${chromosome}_${highlightMin0}_${highlightMax1}.txt";
open (my $out_matrix, '>', $fileMatrix) or die "$!";
my $matrixHeader = "ID\t" . join( "\t", @refBasesAT) . "\n";
print $out_matrix $matrixHeader;

# Original m6A BED12 lines, just sorted (should I redo them with the header?)
my $fileSorted = "$outputFolder/included_m6A_bed_tracks_sorted_by_methylation.txt";
open (my $out_sorted, '>', $fileSorted) or die "$!";

# Same order for outputs
my @methsorted = sort { $b->{meanmeth} <=> $a->{meanmeth} } @hashes;
for my $hashref (@methsorted) {
    my $line = $hashref->{line};
    print $out_sorted "$line\n";
    #
    my $name = $hashref->{name};
    my $aref = $hashref->{vector};
    #print $out_matrix join("\t", ("ZMW$name", @$aref)) . "\n";
    # The name begins with the "m*" movie name so "ZMW" prefix not needed
    print $out_matrix join("\t", ("$name", @$aref)) . "\n";
}

#
close $out_sorted;
close $out_matrix;


exit;
# This mean methylation distribution was only confusing
#my $fileDist = "$outputFolder/mean_methylation.txt";
#open (my $out_dist, '>', $fileDist) or die "$!";
# (loop)
    #my $fraction = $hashref->{meanmeth};
    #print $out_dist "$fraction\n";
#close $out_dist;

