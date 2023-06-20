package chipseq;
use strict;
use Bio::SeqIO;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);
use lib "$Bin";
#use MutationTools::color_space;
use MutationTools::read_Tab;
use MutationTools::read_FASTA;
use MutationTools::read_SAM;
use Data::Dumper;
use lib "/home/unix/rfarrer/perl5/lib/perl5/x86_64-linux-thread-multi/";
use Bio::DB::HTS;

### rfarrer@broadinstitute.org

# try this instead of the bam one to try and check and get consistant read counts with sam
sub count_reads_in_sam_over_peaks {
	my ($file, $peaks, $verbose) = @_;

	my ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak, $unaligned) = (0, 0, 0);


	# Parse SAM file 2
	warn "count_reads_in_sam_over_peaks: Parsing $file\n";
	open my $fh, '<', $file or die "Cannot open $file $!";
	SAM: while (my $line=<$fh>) {
		chomp $line;
		my ($SAM_line) = samlines::read_SAM_lines($line);	

		# Ignore headers
		next SAM if($$SAM_line{'next'} eq 1);
		
		# unaligned
		if($$SAM_line{'aligned'} eq 'N') {
			$unaligned++;
			next SAM;
		}

		my $contig = $$SAM_line{'Ref_seq_name'};
		my $left_pos = $$SAM_line{'leftmost_position'};
		my $right_pos = $$SAM_line{'rightmost_position'};
		
		# Count overlap of peaks
		my %start_to_end_aln;
		$start_to_end_aln{$left_pos}{$right_pos} = 1;
		
		# Count reads that align to peak regions on this contig and count reads that do not
		my ($number_reads_in_region, $number_reads_not_in_region) = &count_reads_in_region($$peaks{$contig}, \%start_to_end_aln);
		$sum_reads_align_in_peak += $number_reads_in_region;
		$sum_reads_align_not_in_peak += $number_reads_not_in_region;

		# tally
		my $total_reads_processed = ($unaligned + $sum_reads_align_in_peak + $sum_reads_align_not_in_peak);
		if($total_reads_processed % 1000000 == 0) {
			warn "count_reads_in_sam_over_peaks: processed $total_reads_processed reads\n";
		}
	}
	close $fh;

	return ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak, $unaligned);
}

# this method seems to suffer from not considering ALL reads in SAM, and also varies slightly with window length (fewer reads for longer window length). 
# Window length isn't really necessary currently, unless this was to be parallelised.
sub count_reads_in_bam_over_peaks {
	my ($file, $fasta, $peaks, $window_length, $verbose) = @_;

	# Save FASTA
	my $sequences = fastafile::fasta_to_struct($fasta);

	# Go through BAM
	my ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak) = (0, 0);
	my $sam = Bio::DB::HTS->new(-bam => $file, -fasta => $fasta);

	# Foreach contig
	CONTIGS: foreach my $contig(sort keys %{$$sequences{'seq'}}) {
		warn "\tProcessing $contig...\n";
		my $contig_length = $$sequences{'seq_length'}{$contig};

		# Foreach non-overlapping window
		my $last_position_in_loop = 1;
		for(my $i=$window_length; $i < $contig_length; $i+=$window_length) {
			my $loop_stop = ($i - 1);

			# For print windows
			my $start_window = $last_position_in_loop;
			if($last_position_in_loop eq 1) { $start_window = ($last_position_in_loop - 1); }

			my $stop_window = ($i + $window_length);
			#if($verbose ne 'n') { warn "Processing $contig $last_position_in_loop - $loop_stop\n"; }
			my @alignments = $sam->get_features_by_location(-seq_id =>$contig,-start => $last_position_in_loop,-end => $loop_stop);
			my ($start_to_end_aln, $total_reads) = &return_alignment_from_BioDBHTS_array(\@alignments);

			if($verbose ne 'n') { warn "Processing $contig $last_position_in_loop - $loop_stop = $total_reads reads\n"; }

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

		if($verbose ne 'n') { warn "Processing $contig $last_position_in_loop - $contig_length = $total_reads reads\n"; }

		# Count reads that align to peak regions on this contig and count reads that do not
		my ($number_reads_in_region, $number_reads_not_in_region) = &count_reads_in_region($$peaks{$contig}, $start_to_end_aln);

		# check
		my $num_reads_counted = ($number_reads_in_region + $number_reads_not_in_region); 
		if($num_reads_counted ne $total_reads) {
			warn "WARNING: Not all reads accounted for (or counted twice) in $contig $last_position_in_loop - $contig_length : reads from BAM parser = $total_reads and num counted = $num_reads_counted\n";
		}

		$sum_reads_align_in_peak += $number_reads_in_region;
		$sum_reads_align_not_in_peak += $number_reads_not_in_region;
	}

	return ($sum_reads_align_in_peak, $sum_reads_align_not_in_peak);
}

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

sub peaks_file_to_chr_start_end_to_one_hash {
	my ($file) = @_;
	warn "save_peaks: $file\n";

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

		#warn "$chr -> $start -> $end = 1\n";
		$peaks{$chr}{$start}{$end} = 1;
	}

	return \%peaks;
}

sub print_peaks_dataframe_for_ggplot {
	my ($peaks_struct, $contig_cumulative_starts, $outfile, $optional_contig, $optional_start, $optional_end, $optional_chipseq_struct) = @_;

	my $max_FC = 0;
	if(-e $outfile) {
		warn "print_peaks_dataframe_for_ggplot: $outfile already made, so not overwriting. Reading it for depth ymax...\n";
		open my $fh, '<', $outfile or die "Error: Cannot open $outfile : $!";
		while(my $line = <$fh>) {
			chomp $line;

			# ignore header
			next if($line =~ m/^file\tpeakID/);

			my @bits = split /\t/, $line;
			my ($file, $tpeakID, $start, $end, $strand, $orientation, $fold_change) = (@bits);

			# max fold change
    		if($fold_change > $max_FC) { $max_FC = $fold_change; }
		}
		return $max_FC;
	}

	warn "print_peaks_dataframe_for_ggplot: printing to $outfile\n";

	# header
	open my $ofh, '>', $outfile or die "Cannot open $outfile : $!";
	print $ofh "file\tpeakID\tstart\tend\tstrand\torientation\tFoldChange\tp_value\n";

	# dataframe
	my $count_of_peaks = 0;
	foreach my $file(sort keys %{$$peaks_struct{'file_to_contig_to_peakID'}}) {
    		foreach my $contig(sort keys %{$$peaks_struct{'file_to_contig_to_peakID'}{$file}}) {

    			# Get cumulative length
    			die "No cumulative length saved for $contig\n" if(!defined $$contig_cumulative_starts{$contig});
    			my $contig_cumulative_start = $$contig_cumulative_starts{$contig};

 	   		# label instead of file
    			my $label = $file;
    			if(defined $$peaks_struct{'peak_file_to_name'}{$file}) { $label = $$peaks_struct{'peak_file_to_name'}{$file}; }

			foreach my $PeakID(keys %{$$peaks_struct{'file_to_contig_to_peakID'}{$file}{$contig}}) {

				$count_of_peaks++;
				my $start       = $$peaks_struct{'file_to_peakID_to_start'}{$file}{$PeakID};
				my $end         = $$peaks_struct{'file_to_peakID_to_end'}{$file}{$PeakID};
				my $strand      = $$peaks_struct{'file_to_peakID_to_strand'}{$file}{$PeakID};
				my $fold_change = $$peaks_struct{'file_to_peakID_to_FC'}{$file}{$PeakID};
				my $pvalue      = $$peaks_struct{'file_to_peakID_to_p_value'}{$file}{$PeakID};

				my $orientation = 1;
				if($strand eq '-') { $orientation = '-1'; }

				if($fold_change > $max_FC) { $max_FC = $fold_change; }

				my $cumulative_start = $start + $contig_cumulative_start;
				my $cumulative_stop = $end + $contig_cumulative_start;

				print $ofh "$label\t$PeakID\t$cumulative_start\t$cumulative_stop\t$strand\t$orientation\t$fold_change\t$pvalue\n";
			}

			# make dummy data if no peaks are given
   			if(($count_of_peaks eq 0) && (defined $optional_start) && (defined $optional_end)) {
   				warn "print_peaks_dataframe_for_ggplot: no peaks found. Creating dummy peak1\n";
   
   				$count_of_peaks++;
				my $mid_point = $contig_cumulative_start + $optional_start + (($optional_end - $optional_start) / 2);
				print $ofh "$label\tdummy_peak1\t$mid_point\t$mid_point\t+\t1\t0\t0\n";
   			}
   		}
   	}

   	# maybe still no peaks (no files at all!)
   	if(($count_of_peaks eq 0) && (defined $optional_start) && (defined $optional_end)) {
   		warn "print_peaks_dataframe_for_ggplot: no peaks found. Creating dummy peak2\n";

   		$count_of_peaks++;
   		die "print_peaks_dataframe_for_ggplot: Error: no peaks found, and no contig specified to make dummy data : $optional_contig" if(!defined $$contig_cumulative_starts{$optional_contig});
   		my $contig_cumulative_start = $$contig_cumulative_starts{$optional_contig};
   		my $mid_point = $contig_cumulative_start + $optional_start + (($optional_end - $optional_start) / 2);

   		#print Dumper($contig_cumulative_starts);
   		#die "end here: $optional_contig, $contig_cumulative_start, $optional_start, $optional_end\n";

		# get label
		my $label;
		NAME:foreach my $labels(keys %{$optional_chipseq_struct}) {
			$label = $labels;
			last NAME;
		}

   		print $ofh "$label\tdummy_peak2\t$mid_point\t$mid_point\t+\t1\t0\t0\n";
   	}
   	

	return $max_FC;
}

sub save_homer_peaks_from_Name_Type_Location {
	my ($Name_Type_Location) = $_[0];

	warn "save_homer_peaks_from_Name_Type_Location...\n";
	my $peaks_struct;

	# save info
	foreach my $name(sort keys %{$Name_Type_Location}) {
		TYPE: foreach my $type(keys %{$$Name_Type_Location{$name}}) {
			
			# only process peaks
			next TYPE if($type ne 'PEAK');

			foreach my $file_location(keys %{$$Name_Type_Location{$name}{$type}}) {
				
				# save metadata
				$$peaks_struct{'number_of_peaks_files'}++;
				if(!defined $$peaks_struct{'first_file'}) { 
					warn "save_homer_peaks_from_Name_Type_Location: first file = $file_location\n";
					$$peaks_struct{'first_file'} = $file_location; 
				}
				$$peaks_struct{'peak_file_to_name'}{$file_location} = $name;

				# save peaks
				$peaks_struct = &save_peaks_either_homer_or_macs2_broad($file_location, $peaks_struct);
			}
		}
	}
	warn "save_homer_peaks_from_Name_Type_Location: Finished\n";
	return ($peaks_struct);
}

sub save_multiple_homer_peaks_files_wrapper {
	my ($homer_files_separated_by_comma) = @_;

	warn "save_multiple_homer_peaks_files_wrapper...\n";
	my $peaks_struct;
	my $first_peak_file;
	my $number_of_peaks_files = 0;
	if($homer_files_separated_by_comma) {
		my @peak_files = split /,/, $homer_files_separated_by_comma;
		$number_of_peaks_files = scalar(@peak_files);

		# Save peaks into memory
		foreach my $peak_file(@peak_files) {
			if(!defined $first_peak_file) { $first_peak_file = $peak_file; }
			$peaks_struct = &save_homer_peaks($peak_file, $peaks_struct);
		}
	}
	return ($peaks_struct, $first_peak_file, $number_of_peaks_files);
}

sub save_peaks_either_homer_or_macs2_broad {
	my ($file, $peaks_struct) = @_;

	my $peak_program = 'macs2';
	warn "save_peaks_either_homer_or_macs2_broad: $file (default assumption $peak_program)\n";

	# Save file
	open my $fh, '<', $file or die "Cannot open $file : $!";
	LINES: while(my $line = <$fh>) {
		chomp $line;

		# ignore header
		if($line =~ m/^#/) {
			if($line =~ m/HOMER Peaks/) {
				warn "save_peaks_either_homer_or_macs2_broad: Homer peaks detected\n";
				$peak_program = 'homer';
			}
			next LINES;
		}

		# check format
		my @bits = split /\t/, $line;

		if($peak_program eq 'homer') {
			die "Error: unexpected Homer peak format: $line\n" if(!defined $bits[12]);
			my ($PeakID, $chr, $start, $end, $strand, $Norm_Tag_Count, $region_size, $findPeaks_score, $Total_Tags, $Control_Tags_norm_to_IP, $Fold_Change_vs_Control, $pvalue_vs_Control, $Clonal_Fold_Change) = @bits;

			# save
			$$peaks_struct{'file_to_contig_to_peakID'}{$file}{$chr}{$PeakID} = 1;
			$$peaks_struct{'file_to_peakID_to_start'}{$file}{$PeakID} = $start;
			$$peaks_struct{'file_to_peakID_to_end'}{$file}{$PeakID} = $end;
			$$peaks_struct{'file_to_peakID_to_strand'}{$file}{$PeakID} = $strand;
			$$peaks_struct{'file_to_peakID_to_FC'}{$file}{$PeakID} = $Fold_Change_vs_Control;
			$$peaks_struct{'file_to_peakID_to_p_value'}{$file}{$PeakID} = $pvalue_vs_Control;
		} else {
			# check format
			# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html
			# and https://nbisweden.github.io/workshop-archive/workshop-ChIP-seq/2018-11-07/labs/lab-broadpeaks.html (broad = narrow minus 10th column)
			my @bits = split /\t/, $line;
			#die "Error: unexpected macs2 broad peak format: $line\n" if((!defined $bits[8]) || (defined $bits[9]));
			die "Error: unexpected macs2 broad peak format: $line\n" if(!defined $bits[8]);
			my ($chr, $start, $end, $PeakID, $score, $strand, $signal_value, $pvalue_vs_Control, $qvalue_vs_Control) = @bits;
			#my ($PeakID, $chr, $start, $end, $strand, $Norm_Tag_Count, $region_size, $findPeaks_score, $Total_Tags, $Control_Tags_norm_to_IP, $Fold_Change_vs_Control, $pvalue_vs_Control, $Clonal_Fold_Change) = @bits;

			# save
			$$peaks_struct{'file_to_contig_to_peakID'}{$file}{$chr}{$PeakID} = 1;
			$$peaks_struct{'file_to_peakID_to_start'}{$file}{$PeakID} = $start;
			$$peaks_struct{'file_to_peakID_to_end'}{$file}{$PeakID} = $end;
			$$peaks_struct{'file_to_peakID_to_strand'}{$file}{$PeakID} = $strand;
			$$peaks_struct{'file_to_peakID_to_FC'}{$file}{$PeakID} = $score; # not quite right. unsure what score is, but probably not FC
			$$peaks_struct{'file_to_peakID_to_p_value'}{$file}{$PeakID} = $pvalue_vs_Control;
		}
	}
	return $peaks_struct;
}

sub peaks_struct_limit_by_contig_start_and_end {
	my ($peaks_struct, $contig_of_interest, $start_of_interest, $end_of_interest) = @_;
	my %new_peaks_struct;

	# transfer metadata
	if(defined $$peaks_struct{'number_of_peaks_files'}) { $new_peaks_struct{'number_of_peaks_files'} = $$peaks_struct{'number_of_peaks_files'}; }
	if(defined $$peaks_struct{'first_file'}) { $new_peaks_struct{'first_file'} = $$peaks_struct{'first_file'}; }

	warn "peaks_struct_limit_by_contig_start_and_end: $contig_of_interest, $start_of_interest, $end_of_interest\n";

	FILE: foreach my $file(sort keys %{$$peaks_struct{'file_to_contig_to_peakID'}}) {

		# transfer metadata
		if(defined $$peaks_struct{'peak_file_to_name'}{$file}) { $new_peaks_struct{'peak_file_to_name'}{$file} = $$peaks_struct{'peak_file_to_name'}{$file}; }

		CONTIG: foreach my $contig(keys %{$$peaks_struct{'file_to_contig_to_peakID'}{$file}}) {
			next if($contig ne $contig_of_interest);

			foreach my $PeakID(keys %{$$peaks_struct{'file_to_contig_to_peakID'}{$file}{$contig}}) {

				my $start_position = $$peaks_struct{'file_to_peakID_to_start'}{$file}{$PeakID};
				my $end_position = $$peaks_struct{'file_to_peakID_to_end'}{$file}{$PeakID};
				my $strand = $$peaks_struct{'file_to_peakID_to_strand'}{$file}{$PeakID};
				my $FC = $$peaks_struct{'file_to_peakID_to_FC'}{$file}{$PeakID};
				my $pvalue = $$peaks_struct{'file_to_peakID_to_p_value'}{$file}{$PeakID};

				if(defined $start_of_interest) {
					next if($start_position < $start_of_interest);
				}
				if(defined $end_of_interest) {
					next if($end_position > $end_of_interest);
				}

				# save
				$new_peaks_struct{'file_to_contig_to_peakID'}{$file}{$contig}{$PeakID} = 1;
				$new_peaks_struct{'file_to_peakID_to_start'}{$file}{$PeakID} = $start_position;
				$new_peaks_struct{'file_to_peakID_to_end'}{$file}{$PeakID} = $end_position;
				$new_peaks_struct{'file_to_peakID_to_strand'}{$file}{$PeakID} = $strand;
				$new_peaks_struct{'file_to_peakID_to_FC'}{$file}{$PeakID} = $FC;
				$new_peaks_struct{'file_to_peakID_to_p_value'}{$file}{$PeakID} = $pvalue;
			}
		}
	}

	return \%new_peaks_struct;
}

1;