#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
#use MutationTools::read_BAM;
use MutationTools::read_FASTA;
use MutationTools::read_Tab;
use MutationTools::read_Chipseq_peaks;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -b <Name tab BAM tab location> -f <fasta reference> -p <Name tab PEAK tab peaks location>\n
Optional: -w Windows for screening genome [10000]
          -v Verbose mode (y/n) [n]
Notes: Does not count unaligned reads
Notes2: print to peak file-FRiP-scores-RF.txt
Notes3: Works for either macs2 or homer\n";
our($opt_b, $opt_f, $opt_p, $opt_v, $opt_w);
getopt('bfpvw');
die $usage unless ($opt_b && $opt_f && $opt_p);
foreach($opt_b, $opt_f, $opt_p) { die "error: Cannot open $_ : $!\n" unless (-e $_); }
if(!defined $opt_v) { $opt_v = 'n'; }
if(!defined $opt_w) { $opt_w = 10000; }

# Save Name to BAM location
my $BAM_locations = tabfile::save_columns_to_column_hash($opt_b, 0, 2);

# Save Name to Peak locations
my $Peak_locations = tabfile::save_columns_to_column_hash($opt_p, 0, 2);

# sanity check
foreach my $name(keys %{$Peak_locations}) {
	die "No $name in both BAM and Peaks\n" if(!defined $$BAM_locations{$name});
}

# go through them.
foreach my $name(keys %{$Peak_locations}) {
	my $bam_file = $$BAM_locations{$name};
	my $peak_file = $$Peak_locations{$name};

	# Save peaks
	my $peaks = chipseq::peaks_file_to_chr_start_end_to_one_hash($opt_p);


	# Go through BAM
	my ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak) = chipseq::count_reads_in_bam_over_peaks($bam_file, $opt_f, $peaks, $opt_w, $opt_v);


	# Summary
	my $sum_reads_found = ($sum_reads_align_in_peak + $sum_reads_align_not_in_peak);
	my $FRIP = ($sum_reads_align_in_peak / $sum_reads_found)*100;

	# Output
	my $out_file = "$peak_file-FRiP-scores-RF.txt";
	open my $ofh, '>', $out_file or die "Cannot open $out_file : $!";
	print $ofh "command_run:$0 -b $opt_b -f $opt_f -p $opt_p\n";
	print $ofh "bam processed:\t$bam_file\n";
	print $ofh "peaks processed:'t'$peak_file\n";
	print $ofh "Reads in peaks\t$sum_reads_align_in_peak\n";
	print $ofh "Reads not in peaks\t$sum_reads_align_not_in_peak\n";
	print $ofh "FRiP score\t$FRIP\n";
}