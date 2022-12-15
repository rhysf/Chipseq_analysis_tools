#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
#use MutationTools::read_BAM;
use read_FASTA;
use Data::Dumper;
use lib "/home/unix/rfarrer/perl5/lib/perl5/x86_64-linux-thread-multi/";
use Bio::DB::HTS;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -b <BAM> -f <fasta reference> -p <peaks file>\n
Optional: -t Peak caller (homer/macs2) [homer]
          -w Windows for screening genome [10000]
Notes: Does not count unaligned reads\n";
our($opt_b, $opt_f, $opt_p, $opt_t, $opt_w);
getopt('bfptw');
die $usage unless ($opt_b && $opt_f && $opt_p);
foreach($opt_b, $opt_f, $opt_p) { die "error: Cannot open $_ : $!\n" unless (-e $_); }
if(!defined $opt_t) { $opt_t = 'homer'; }
if(!defined $opt_w) { $opt_w = 10000; }

# Save peaks
my $peaks = &save_peaks($opt_p, $opt_t);

# Save FASTA
my $sequences = fastafile::fasta_to_struct($opt_f);

# Go through BAM
my $sam = Bio::DB::HTS->new(-bam => $opt_b, -fasta=> $opt_f);

my ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak) = (0, 0);


# Foreach contig
CONTIGS: foreach my $contig(sort keys %{$$sequences{'seq_length'}}) {
	warn "\tProcessing $contig...\n";
	my $contig_length = $$sequences{'seq_length'}{$contig};

	# Foreach non-overlapping window
	my $last_position_in_loop = 1;
	for(my $i=$opt_w; $i < $contig_length; $i+=$opt_w) {
		my $loop_stop = ($i - 1);

		# For print windows
		my $start_window = $last_position_in_loop;
		if($last_position_in_loop eq 1) { $start_window = ($last_position_in_loop - 1); }

		my $stop_window = ($i + $opt_w);
		#warn "Processing $contig $last_position_in_loop - $loop_stop\n";
		my @alignments = $sam->get_features_by_location(-seq_id =>$contig,-start => $last_position_in_loop,-end => $loop_stop);
		my ($start_to_end_aln, $total_reads) = &return_alignment_from_BioDBHTS_array(\@alignments);

		#print Dumper($start_to_end_aln);
		#die "end here\n";

		# Count reads that align to peak regions on this contig and count reads that do not
		my ($number_reads_in_region, $number_reads_not_in_region) = &count_reads_in_region($$peaks{$contig}, $start_to_end_aln);

		# check
		my $num_reads_counted = ($number_reads_in_region + $number_reads_not_in_region); 
		if($num_reads_counted ne $total_reads) {
			warn "WARNING: Not all reads accounted for (or counted twice) in $contig $last_position_in_loop - $loop_stop : reads from BAM parser = $total_reads and num counted = $num_reads_counted\n";
		}

		$sum_reads_align_in_peak += $number_reads_in_region;
		$sum_reads_align_not_in_peak += $number_reads_not_in_region;

		#print "$contig from $last_position_in_loop to $loop_stop has $number_reads_in_region reads in peaks and $number_reads_not_in_region reads not in peaks\n";

		$last_position_in_loop = $i;
	}

	# Last group of reads for contig
	#warn "Processing $contig $last_position_in_loop - $$fasta_lengths{$contig}\n";
	my $start_window = $last_position_in_loop;
	if($last_position_in_loop eq 1) { $start_window = ($last_position_in_loop - 1); }
	my @alignments = $sam->get_features_by_location(-seq_id =>$contig,-start => $last_position_in_loop,-end => $contig_length);
	my ($start_to_end_aln, $total_reads) = &return_alignment_from_BioDBHTS_array(\@alignments);

	# Count reads that align to peak regions on this contig and count reads that do not
	my ($number_reads_in_region, $number_reads_not_in_region) = &count_reads_in_region($$peaks{$contig}, $start_to_end_aln);
	$sum_reads_align_in_peak += $number_reads_in_region;
	$sum_reads_align_not_in_peak += $number_reads_not_in_region;
}

#output
my $sum_reads_found = ($sum_reads_align_in_peak + $sum_reads_align_not_in_peak);
my $FRIP = ($sum_reads_align_in_peak / $sum_reads_found)*100;
print "$0 -b $opt_b -f $opt_f -p $opt_p\n";
print "Reads in peaks\t$sum_reads_align_in_peak\n";
print "Reads not in peaks\t$sum_reads_align_not_in_peak\n";
print "FRiP score\t$FRIP\n";

sub count_reads_in_region {
	my ($peaks, $aln) = @_;

	#print Dumper($peaks);
	#die;

	my ($reads_in_peaks, $reads_not_in_peaks) = (0, 0);
	my $reads_ive_processed = 0;

	# Go through each alignment
	foreach my $start(keys %{$aln}) {
		foreach my $end(keys %{$$aln{$start}}) {
			my $number_reads = $$aln{$start}{$end};
			$reads_ive_processed += $number_reads;

			#warn "--->read $start - $end - $number_reads\n";

			# Found in a peak?
			my $found = 0;

			# Go through each peak
			PEAK_START: foreach my $peak_start(keys %{$peaks}) {
				PEAK_END: foreach my $peak_end(keys %{$$peaks{$peak_start}}) {

					#warn "checking peak start $peak_start and peak end $peak_end\n";
					
					LOOP: for(my $i=$start; $i<$end; $i++) {
						if(($i >= $peak_start) && ($i <= $peak_end)) {
							#warn "\treads in peaks += $number_reads\n";
							$reads_in_peaks += $number_reads;
							$found = 1;
							last PEAK_START;
						}
					}

					#warn "does read at $start - $end ($number_reads) fall within peak $peak_start - $peak_end\n";
					#die;
				}
			}

			# Not found in a peak
			if($found eq 0) { 
				#warn "\treads not in peaks += $number_reads\n";
				$reads_not_in_peaks += $number_reads; 
			}
		}
	}
	#warn "count_reads_in_region: processed $reads_ive_processed reads = $reads_in_peaks reads in peaks + $reads_not_in_peaks reads not in peaks\n";
	return ($reads_in_peaks, $reads_not_in_peaks);
}

sub save_peaks {
	my ($file, $peak_caller) = @_;
	warn "save_peaks: $file, $peak_caller\n";

	my %peaks;
	open my $fh, '<', $file or die "error: Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		next if($line =~ m/^#/);
		my @bits = split /\t/, $line;

		my ($chr, $start, $end);

		# Peak caller
		# homer
		if($peak_caller eq 'homer') {
			($chr, $start, $end) = ($bits[1], $bits[2], $bits[3]);
		}
		# macs2
		elsif($peak_caller eq 'macs2') {
			($chr, $start, $end) = ($bits[0], $bits[1], $bits[2]);
		}
		else {
			die "error: unrecognised peak caller: $peak_caller\n";
		}

		#warn "$chr -> $start -> $end = 1\n";
		$peaks{$chr}{$start}{$end} = 1;
	}

	return \%peaks;
}

sub return_alignment_from_BioDBHTS_array {
	my ($alignments) = $_[0];
	my %saved_regions;
	my $total_reads_saved = 0;
	SAMSEQS: for my $a(@{$alignments}) {
		# where does the alignment start in the reference sequence
		my $start = $a->start;
		my $end = $a->end;
		# query/read sequence
		my $query_dna = $a->query->dna; 
		# where does the alignment start in the query sequence
		my $query_start = $a->query->start;     
		my $query_end   = $a->query->end;
		# new for determining indels
		my $cigar = $a->cigar_str;

		# Ignore reads specifying/aligning over indels
		#my $read_length = ($query_end - $query_start);
		#my $ref_seq_length = ($end - $start);
		#next SAMSEQS if ($read_length ne $ref_seq_length);

		# Save
		#my $line = join "\t", $start, $end, $query_dna, $query_start, $query_end, $cigar;
		#$lines .= "$line\n";
		#warn "query start $query_start - query end $query_end and start $start and end $end\n";
		$saved_regions{$start}{$end}++;
		$total_reads_saved++;
	}
	return (\%saved_regions, $total_reads_saved);
}
