#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
#use FindBin qw($Bin);
#use lib "$Bin/perl_modules";

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -a <peak calls1> -b <peak calls 2>\n";
our($opt_a, $opt_b);
getopt('ab');
die $usage unless ($opt_a && $opt_b);
foreach my $file($opt_a, $opt_b) { die "Cannot open $file : $!" unless (-e $file); }

# save peak calls (different peaks can be given the same name across multiple file e.g. 7000000096786010-222 in file 1 and file2 - so give a unique tag)
my $peak_calls1 = &peaks_file_to_struct($opt_a);
my $peak_calls2 = &peaks_file_to_struct($opt_b);

# compare bidirectional
my $overlap1_and_2 = &compare($peak_calls1, $peak_calls2);

# print header
print "Peak1\tPeak2\n";

#print Dumper($overlap1_and_2);
#die;

# print file 1 with overlaps
foreach my $peak1(sort keys %{$peak_calls1}) {
	my $peak1_line = $$peak_calls1{$peak1}{'line'};

	# overlaps saved?
	my $overlap_entry1;
	if(defined $$overlap1_and_2{$peak1}) {
		foreach my $peaks_in_file2(keys %{$$overlap1_and_2{$peak1}}) {
			$overlap_entry1 .= ($peaks_in_file2 . ', ');
		}
		$overlap_entry1 =~ s/,\s$//;
	}

	if(!defined $overlap_entry1) { $overlap_entry1 = "None"; }
	else {

		print "$peak1\t$overlap_entry1\n";
	}
}

sub compare {
	my ($info1, $info2) = @_;
	my %overlapping_entries;

	# go through every peak in file 1
	warn "compare: files 1 to 2...\n";
	PEAK1: foreach my $peak1(keys %{$info1}) {
		my $peak1_chr = $$info1{$peak1}{'chr'};
		my $peak1_start = $$info1{$peak1}{'start'};
		my $peak1_end = $$info1{$peak1}{'end'};

		# overlap with any in file2?
		PEAK2: foreach my $peak2(keys %{$info2}) {
			my $peak2_chr = $$info2{$peak2}{'chr'};
			my $peak2_start = $$info2{$peak2}{'start'};
			my $peak2_end = $$info2{$peak2}{'end'};

			# ignore if not on same chromosome/contig
			next if($peak1_chr ne $peak2_chr);

			# check if overlap? (slow)
			FIND: for(my $i=$peak1_start; $i<=$peak1_end; $i++) {
				if(($i eq $peak2_start) || ($i eq $peak2_end)) {
					#warn "1-2>: saving $peak1 -> $peak2 (and vice versa)\n";
					$overlapping_entries{$peak1}{$peak2} = 1;
					$overlapping_entries{$peak2}{$peak1} = 1;
					#last PEAK2;
				}
			}
		}
	}

	warn "compare: files 2 to 1...\n";
	PEAK2: foreach my $peak2(keys %{$info2}) {
		my $peak2_chr = $$info2{$peak2}{'chr'};
		my $peak2_start = $$info2{$peak2}{'start'};
		my $peak2_end = $$info2{$peak2}{'end'};

		# overlap with any in file2?
		PEAK1: foreach my $peak1(keys %{$info1}) {
			my $peak1_chr = $$info1{$peak1}{'chr'};
			my $peak1_start = $$info1{$peak1}{'start'};
			my $peak1_end = $$info1{$peak1}{'end'};

			# ignore if not on same chromosome/contig
			next if($peak1_chr ne $peak2_chr);

			# ignore if already saved
			#next if(exists $overlapping_entries{$peak1}{$peak2});
			#next if(exists $overlapping_entries{$peak2}{$peak1});

			#warn "A new unaccounted peak was found from file 2 -> 1: $peak2 -> $peak1\n";

			# check if overlap? (slow)
			FIND: for(my $i=$peak1_start; $i<=$peak1_end; $i++) {
				if(($i eq $peak2_start) || ($i eq $peak2_end)) {
					#warn "2-1>: saving $peak1 -> $peak2 (and vice versa)\n";
					$overlapping_entries{$peak1}{$peak2} = 1;
					$overlapping_entries{$peak2}{$peak1} = 1;
					#last PEAK2;
				}
			}
		}
	}

	#die "end here\n";

	return \%overlapping_entries;
}

sub peaks_file_to_struct {
	my ($file) = @_;
	warn "peaks_file_to_struct: $file\n";

	my $peak_caller = "macs2";

	my %peaks;
	open my $fh, '<', $file or die "error: Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;

		# ignore headers
		if($line =~ m/# HOMER Peaks/) {
			warn "save_peaks: Homer peak caller detected.\n";
			$peak_caller = 'homer';
		}
		next if($line =~ m/^#/);
		next if($line =~ m/^\n/);
		next if($line eq '');

		my @bits = split /\t/, $line;

		my ($chr, $start, $end);

		# Peak callers
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

		# save (peak name -> chr, start, end, fold_change and pvalue)
		my $PeakID = ($chr . '-' . $start . '-' . $tag);
		$info{$PeakID}{'chr'} = $chr;
		$info{$PeakID}{'start'} = $start;
		$info{$PeakID}{'end'} = $end;
		$info{$PeakID}{'line'} = $line;
	}

	return \%peaks;
}