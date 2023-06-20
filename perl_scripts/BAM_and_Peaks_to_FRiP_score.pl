#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
#use MutationTools::read_BAM;
use MutationTools::read_Chipseq_peaks;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -b <BAM> -f <fasta reference> -p <peaks file>\n
Optional: -w Windows for screening genome [10000]
          -v Verbose mode (y/n) [n]
Notes: Does not count unaligned reads
Notes2: Works for either macs2 or homer\n";
our($opt_b, $opt_f, $opt_p, $opt_v, $opt_w);
getopt('bfptvw');
die $usage unless ($opt_b && $opt_f && $opt_p);
foreach($opt_b, $opt_f, $opt_p) { die "error: Cannot open $_ : $!\n" unless (-e $_); }
if(!defined $opt_v) { $opt_v = 'n'; }
if(!defined $opt_w) { $opt_w = 10000; }

# Save peaks
my $peaks = chipseq::peaks_file_to_chr_start_end_to_one_hash($opt_p);

# Go through BAM
my ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak) = chipseq::count_reads_in_bam_over_peaks($opt_b, $opt_f, $peaks, $opt_w, $opt_v);

#output
my $sum_reads_found = ($sum_reads_align_in_peak + $sum_reads_align_not_in_peak);
my $FRIP = ($sum_reads_align_in_peak / $sum_reads_found)*100;
print "$0 -b $opt_b -f $opt_f -p $opt_p\n";
print "Reads in peaks\t$sum_reads_align_in_peak\n";
print "Reads not in peaks\t$sum_reads_align_not_in_peak\n";
print "FRiP score\t$FRIP\n";