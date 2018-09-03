#!/usr/bin/env perl

# comp_freq.pl is Perl script for computing codon count, amino acid 
# count, codon pair count, amino acid pair count, codon frequency 
# (ICU), and codon context (CC) from a fasta file.
#
# For the definition of individual codon usage (ICU) and codon context 
# (CC), please see: Chung BK, Lee DY, Computational codon optimization 
# of synthetic gene for protein expression. BMC Syst Biol, 2012, 6:134.
#
# Accepted input is a fasta formated file containing nucleotide 
# sequences (A/C/G/T/U) without any ambiguous bases.
# Computes the codon count, amino acid count, codon pair count, amino 
# acid pair count, codon frequency (ICU), and codon context (CC).
#

use strict;
use warnings;

my %c2aa = ();
my %aa2c = ();
my %codon_count = ();
my %codonp_count = ();
my %cfreq_count = ();
my %cpfreq_count = ();
my %aa_count = ();
my %aap_count = ();

sub read_trans_table
{
	my $fileh;

	open $fileh, "<"."trans_table.txt" or die $!;
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
}

sub init_hashes
{
	foreach my $codon (keys %c2aa) {
		$codon_count{$codon} = 0;
		$cfreq_count{$codon} = 0.0;
		if ($c2aa{$codon} eq "*") {
			next;
		}
		foreach my $codon2 (keys %c2aa) {
			$codonp_count{$codon.$codon2} = 0;
			$cpfreq_count{$codon.$codon2} = 0.0;
		}
	}
	foreach my $aa (keys %aa2c) {
		$aa_count{$aa} = 0;
		if ($aa eq "*") {
			next;
		}
		foreach my $aa2 (keys %aa2c) {
			$aap_count{$aa.$aa2} = 0;
		}
	}
}

sub read_fasta
{
	my $in_filename = shift;

	my $count = 0;
	my $header = "";
	my $seq = "";
	my $fileh;

	open $fileh, "<".$in_filename or die $!;
	my $line = <$fileh>;
	while ($line) {
		$line =~ s/\n|\r//g;
		if ($line =~ /^>/) {
			if ($count) {
				process_seq($header, $seq);
			}
			$header = substr $line, 1;
			$seq = "";
			++$count;
			if (!($count%100)) {
				printf "%d sequences processed.\n", $count;
			}
		} else {
			$line = uc $line;
			$line =~ s/U/T/g;
			$line =~ s/\s//g;
			$seq .= $line;
		}
		$line = <$fileh>;
	}
	process_seq($header, $seq);
	close $fileh;
	printf "Total of %d sequences processed.\n", $count;
}

sub process_seq
{
	my $header = shift;
	my $seq = shift;

	my $internal_stop = 0;
	my $unknown = 0;

	if (length($seq)%3) {
		printf "Warning: length of %s is not divisible by 3.\n", $header;
		$seq = substr($seq, 0, -(length($seq)%3));
	}
	if (length($seq) < 6) {
		printf "Warning: %s is too short to be considered. Skipping.\n", $header;
		return;
	}
	my @stops = @{ $aa2c{"*"} };
	my $rgx = join "|", @stops;
	"a" =~ /a/;
	if ($seq =~ /((?:$rgx){2,})$/) {
		$seq = substr($seq, 0, 3-length($1));
		printf "%s has more than one stop signal at the end.\n", $header;
		
		my @stop_array = unpack "(a3)*", $1;
		foreach my $stopcodon (@stop_array){
			++$codon_count{$stopcodon};
			++$aa_count{$c2aa{$stopcodon}};
		}
	}

	my @codon_array = unpack "(a3)*", $seq;

	if (!exists $c2aa{$codon_array[0]}) {
		++$unknown;
	} elsif ($c2aa{$codon_array[0]} eq "*") {
		++$internal_stop;
	} elsif ($c2aa{$codon_array[0]} ne "M") {
		printf "%s does not start with a methionine.\n", $header;
	} else {
		++$codon_count{$codon_array[0]};
		++$aa_count{$c2aa{$codon_array[0]}};
	}
	for (my $i=0; $i<$#codon_array-1; ++$i) {
		(my $flag1, my $flag2) = process_codons($codon_array[$i], $codon_array[$i+1], 0, $i);
		$internal_stop += $flag1;
		$unknown += $flag2;
	}
	(my $flag1, my $flag2) = process_codons($codon_array[$#codon_array-1], $codon_array[$#codon_array], 1, $#codon_array-1);
	$internal_stop += $flag1;
	$unknown += $flag2;

	if ($internal_stop) {
		printf "%s has internal stop codons.\n", $header;
	}
	if ($unknown) {
		printf "%s has unrecognized codons.\n", $header;
	}
}

sub process_codons
{
	my $codon1 = shift;
	my $codon2 = shift;
	my $second_stop = shift;
	my $pos = shift;

	my $internal_stop = 0;
	my $unknown = 0;

	if (exists $c2aa{$codon2}) {
		if ($c2aa{$codon2} ne "*" || $second_stop) {
			++$codon_count{$codon2};
			++$aa_count{$c2aa{$codon2}};
			if (exists $c2aa{$codon1} && $c2aa{$codon1} ne "*") {
				if ($pos || $c2aa{$codon1} eq "M") {
					++$codonp_count{$codon1.$codon2};
					++$aap_count{$c2aa{$codon1}.$c2aa{$codon2}};
				}
			}
		} else {
			++$internal_stop;
		}
	} else {
		++$unknown;
	}
	return $internal_stop, $unknown;
}

sub print2file
{
	my $fh;

	open $fh, ">"."count_codon.txt" or die $!;
	foreach my $codon (sort{$a cmp $b} keys %c2aa) {
		printf $fh "%s\t%d\n", $codon, $codon_count{$codon};
	}
	close $fh;
	open $fh, ">"."count_codonp.txt" or die $!;
	foreach my $codon1 (sort{$a cmp $b} keys %c2aa) {
		if ($c2aa{$codon1} eq "*") {
			next;
		}
		foreach my $codon2 (sort{$a cmp $b} keys %c2aa) {
			printf $fh "%s\t%d\n", $codon1.$codon2, $codonp_count{$codon1.$codon2};
		}
	}
	close $fh;

	open $fh, ">"."count_aa.txt" or die $!;
	foreach my $aa (sort{$a cmp $b} keys %aa2c) {
		printf $fh "%s\t%d\n", $aa, $aa_count{$aa};
	}
	close $fh;
	open $fh, ">"."count_aap.txt" or die $!;
	foreach my $aa1 (sort{$a cmp $b} keys %aa2c) {
		if ($aa1 eq "*") {
			next;
		}
		foreach my $aa2 (sort{$a cmp $b} keys %aa2c) {
			printf $fh "%s\t%d\n", $aa1.$aa2, $aap_count{$aa1.$aa2};
		}
	}
	close $fh;

	open $fh, ">"."count_icu.txt" or die $!;
	foreach my $codon (sort{$a cmp $b} keys %c2aa) {
		my $aa = $c2aa{$codon};
		if ($codon_count{$codon}) {
			$cfreq_count{$codon} =  $codon_count{$codon} / $aa_count{$aa};
		} else {
			$cfreq_count{$codon} = 0.0;
		}
		printf $fh "%s\t%f\n", $codon, $cfreq_count{$codon};
	}
	close $fh;
	open $fh, ">"."count_cc.txt" or die $!;
	foreach my $codon1 (sort{$a cmp $b} keys %c2aa) {
		if ($c2aa{$codon1} eq "*") {
			next;
		}
		foreach my $codon2 (sort{$a cmp $b} keys %c2aa) {
			my $codonp = $codon1.$codon2;
			my $aap = $c2aa{$codon1}.$c2aa{$codon2};
			if ($codonp_count{$codonp}) {
				$cpfreq_count{$codonp} =  $codonp_count{$codonp} / $aap_count{$aap};
			} else {
				$cpfreq_count{$codonp} = 0.0;
			}
			printf $fh "%s\t%f\n", $codonp, $cpfreq_count{$codonp};
		}
	}
	close $fh;
}

sub main
{
	if ($#ARGV < 0) {
		printf "Usage: comp_freq.pl inputfilename";
		return;
	}
	my $in_filename = $ARGV[0];

	read_trans_table();
	init_hashes();
	read_fasta($in_filename);
	print2file();
}

main();














1
