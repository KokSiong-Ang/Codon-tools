#!/usr/bin/env perl

#
# comp_score.pl is a Perl script for computing ICU score, CC score, 
# and CAI score of sequences with respect to a reference host. It 
# also computes the GC content, GC content at third nucleotide 
# position (GC3), hidden stop count, exclusion sequence count, and 
# consecutive repeats count. The sequences should be supplied in a 
# fasta file. 
#
# For the definition of individual codon usage (ICU) and codon context 
# (CC), please see: Chung BK, Lee DY, Computational codon optimization 
# of synthetic gene for protein expression. BMC Syst Biol, 2012, 6:134.
#
# Accepted input is a fasta formated file containing nucleotide 
# sequences (A/C/G/T/U) without any ambiguous bases.
# Additional input files: 
# Requires trans_table.txt for host codon translation table (see file 
# in folder for format)
# Requires count_icu.txt for host ICU counts (see file in folder for 
# format)
# Requires count_cc.txt for host CC counts (see file in folder for 
# format)
# (You can use the count_icu and count_cc files from comp_freq.pl)
# Optional: exclusion_seq.txt for exclusion sequences (see file in 
# folder for format)
# Optional: repeat_numbers.txt for repeat sequence specification (see 
# file in folder for format)
#
	
use strict;
use warnings;

require "./codon_score_func.pl";


sub read_fasta
{
	my $in_filename = shift;
	my $out_filename = shift;

	my $c2aa_ref = shift;
	my $aa2c_ref = shift;
	my $num_ic = shift;
	my $num_cc = shift;
	my $icu_ref = shift;
	my $cc_ref = shift;
	my $most_freq_ic_value_ref = shift;
	my $exclu_ref = shift;
	my $repeat_ref = shift;

	my $header = "";
	my $seq = "";

	open my $fileh_in, "<".$in_filename or die $!;
	open my $fileh_out, ">".$out_filename or die $!;
	printf $fileh_out "ICU score\tCC score\tCAI score\tHidden\tGC content\tExcluseq\tRepeat\n";

	my $line = <$fileh_in>;
	while ($line) {
		$line =~ s/\n|\r//g;
		if ($line =~ /^>/) {
			if ($seq) {
				process_seq($fileh_out, $header, $seq, $c2aa_ref, $aa2c_ref, $num_ic, $num_cc, $icu_ref, $cc_ref, $most_freq_ic_value_ref, $exclu_ref, $repeat_ref);
			}
			$header = substr $line, 1;
			$seq = "";
		} else {
			$line = uc $line;
			$line =~ s/U/T/g;
			$line =~ s/\s//g;
			$seq .= $line;
		}
		$line = <$fileh_in>;
	}
	process_seq($fileh_out, $header, $seq, $c2aa_ref, $aa2c_ref, $num_ic, $num_cc, $icu_ref, $cc_ref, $most_freq_ic_value_ref, $exclu_ref, $repeat_ref);
	close $fileh_in;
	close $fileh_out;
}

sub process_seq
{
	my $fileh_out = shift;
	my $header = shift;
	my $seq = shift;

	my $c2aa_ref = shift;
	my $aa2c_ref = shift;
	my $num_ic = shift;
	my $num_cc = shift;
	my $icu_ref = shift;
	my $cc_ref = shift;
	my $most_freq_ic_value_ref = shift;
	my $exclu_ref = shift;
	my $repeat_ref = shift;

	if (length($seq)%3) {
		printf "Warning: length of %s is not divisible by 3.\n", $header;
		$seq = substr($seq, 0, -(length($seq)%3));
	}
	if (length($seq) < 6) {
		printf "Warning: %s is too short to be considered. Skipping.\n", $header;
		return;
	}
	my @stops = @{ $aa2c_ref->{"*"} };
	my $rgx = join "|", @stops;
	"a" =~ /a/;
	if ($seq =~ /((?:$rgx){2,})$/) {
		$seq = substr($seq, 0, 3-length($1));
		printf "%s has more than one stop signal at the end.\n", $header;
	}

	(my $c_counts_ref, my $aa_counts_ref) = comp_ic_counts($seq, $header, $c2aa_ref, $aa2c_ref);
	(my $cc_counts_ref, my $aap_counts_ref) = comp_cc_counts($seq, $header, $c2aa_ref, $aa2c_ref);

	my $score_icu = - comp_icu($c2aa_ref, $icu_ref, $num_ic, $c_counts_ref, $aa_counts_ref);
	my $score_cc = - comp_cc($c2aa_ref, $cc_ref, $num_cc, $cc_counts_ref, $aap_counts_ref);
	my $score_cai = - comp_cai($seq, $c2aa_ref, $icu_ref, $most_freq_ic_value_ref);
	my $score_hidden = comp_hidden($seq, $aa2c_ref);

	my $gc_content = comp_gc($seq);
	my $gc3_content = comp_gc3($seq);
	my $score_exc = comp_exc($seq, $exclu_ref);
	my $score_repeat = comp_repeat($seq, $repeat_ref);

	printf $fileh_out "%f\t%f\t%f\t%6d\t%f\t%f\t%6d\t%6d\n", 
			$score_icu, $score_cc, $score_cai, $score_hidden, $gc_content, $gc3_content, 
			$score_exc, $score_repeat;
}

sub read_trans_table
{
	my $trans_table_filename = shift;
	my %c2aa = ();
	my %aa2c = ();

	open my $fileh, "<".$trans_table_filename or die $!;
	my $line = <$fileh>;
	while ($line) {
		$line =~ s/\n|\r//g;
		my @array = split /\t/, $line;
		$array[0] = uc $array[0];
		$array[0] =~ s/U/T/g;
		$array[1] = uc $array[1];
		$array[1] =~ s/\./\*/g;
		$c2aa{$array[0]} = $array[1];
		push @{ $aa2c{$array[1]} }, $array[0];
		$line = <$fileh>;
	}
	close $fileh;
	return \%c2aa, \%aa2c;
}

sub read_icu_table
{
	my $ic_table_filename = shift;
	my $c2aa_ref = shift;
	my %icu_freq = ();

	foreach my $codon (keys %{$c2aa_ref}) {
		$icu_freq{$codon} = 0.0;
	}

	open my $fileh, "<".$ic_table_filename or die $!;
	my $line = <$fileh>;
	while ($line) {
		$line =~ s/\n|\r//g;
		my @array = split /\t/, $line;
		$array[0] = uc $array[0];
		$array[0] =~ s/U/T/g;
		$icu_freq{$array[0]} = $array[1];
		$line = <$fileh>;
	}
	close $fileh;
	return \%icu_freq, scalar keys %icu_freq;
}

sub read_cc_table
{
	my $cc_table_filename = shift;
	my $c2aa_ref = shift;
	my %cc_freq = ();

	foreach my $codon1 (keys %{$c2aa_ref}) {
		if ($c2aa_ref->{$codon1} eq "*") {
			next;
		}
		foreach my $codon2 (keys %{$c2aa_ref}) {
			$cc_freq{$codon1.$codon2} = 0.0;
		}
	}

	open my $fileh, "<".$cc_table_filename or die $!;
	my $line = <$fileh>;
	while ($line) {
		$line =~ s/\n|\r//g;
		my @array = split /\t/, $line;
		$array[0] = uc $array[0];
		$array[0] =~ s/U/T/g;
		$cc_freq{$array[0]} = $array[1];
		$line = <$fileh>;
	}
	close $fileh;
	return \%cc_freq, scalar keys %cc_freq;
}

sub read_exclusion
{
	my $exclu_filename = shift;

	my %exclu = ();
	if (open my $fileh, "<".$exclu_filename) {
		my $line = <$fileh>;
		while ($line) {
			$line =~ s/\n|\r//g;
			$line =~ s/\s//g;
			if ($line =~ /^\s*#/ || length($line) < 1) {
				$line = <$fileh>;
				next;
			}
			$line = uc $line;
			$line =~ s/U/T/g;
			my @array = split ";", $line;
			foreach my $exc (@array) {
				if ($exc !~ /^[ACGT]+$/) {
					printf "Error: non ACGTU characters found in supplied exclusion sequences (%s).\n", $exclu_filename;
				}
				$exclu{$exc} = 1;
			}
			$line = <$fileh>;
		}
		close $fileh;
	} else {
		printf "No %s file found. Assuming no exclusion sequences.\n", $exclu_filename;
	}
	my @exclu_array = keys %exclu;
	return \@exclu_array;
}

sub read_repeat
{
	my $repeat_filename = shift;

	my %repeatnums = ();
	if (open my $fileh, "<".$repeat_filename) {
		my $line = <$fileh>;
		while ($line) {
			$line =~ s/\n|\r//g;
			$line =~ s/\s//g;
			if ($line =~ /^\s*#/ || length($line) < 1) {
				$line = <$fileh>;
				next;
			}
			my @array = split ";", $line;
			foreach my $exc (@array) {
				my @array2 = split ":", $exc;
				$repeatnums{$array2[0]} = $array2[1];
			}
			$line = <$fileh>;
		}
		close $fileh;
	} else {
		printf "No %s file found. Assuming no repeats specified.\n", $repeat_filename;
	}
	return \%repeatnums;
}

sub main
{
	if ($#ARGV < 0) {
		printf "Usage: codon_compute_score.pl inputfilename";
		return;
	}
	my $in_filename = $ARGV[0];
	my $out_filename = "scores.txt";
	my $trans_table_filename = "trans_table.txt";
	my $ic_table_filename = "count_icu.txt";
	my $cc_table_filename = "count_cc.txt";
	my $exclu_filename = "exclusion_seq.txt";
	my $repeat_filename = "repeat_numbers.txt";

	(my $c2aa_ref, my $aa2c_ref) = read_trans_table($trans_table_filename);
	(my $icu_ref, my $num_ic) = read_icu_table($ic_table_filename, $c2aa_ref);
	(my $cc_ref, my $num_cc) = read_cc_table($cc_table_filename, $c2aa_ref);
	my $exclu_ref = read_exclusion($exclu_filename);
	my $repeat_ref = read_repeat($repeat_filename);
	(my $most_freq_ic_ref, my $most_freq_ic_value_ref) = find_most_freq_ic($icu_ref, $c2aa_ref);

	read_fasta($in_filename, $out_filename,	
			$c2aa_ref, $aa2c_ref, 
			$num_ic, $num_cc, 
			$icu_ref, $cc_ref, $most_freq_ic_value_ref, $exclu_ref, $repeat_ref);

}

main();




















