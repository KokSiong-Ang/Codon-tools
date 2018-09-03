#!/usr/bin/env perl

use strict;
use warnings;

sub comp_ic_counts
{
	my $seq = shift;
	my $header = shift;
	my $c2aa_ref = shift;
	my $aa2c_ref = shift;
	my %c_counts = ();
	my %aa_counts = ();

	foreach my $codon (keys %{$c2aa_ref}) {
		$c_counts{$codon} = 0;
	}
	foreach my $aa (keys %{$aa2c_ref}) {
		$aa_counts{$aa} = 0;
	}

	my @codons_array = unpack "(a3)*", $seq;
	my $codon = shift @codons_array;
	if (exists $c2aa_ref->{$codon} && $c2aa_ref->{$codon} eq "M") {
		++$c_counts{$codon};
		++$aa_counts{$c2aa_ref->{$codon}};
	}
	foreach my $codon (@codons_array) {
		if (exists $c2aa_ref->{$codon}) {
			++$c_counts{$codon};
			++$aa_counts{$c2aa_ref->{$codon}};
		} else {
			printf "Unrecognized \"%s\" codon found in %s\n", $codon, $header;
		}
	}
	return \%c_counts, \%aa_counts;
}

sub find_most_freq_ic
{
	my $icu_ref = shift;
	my $c2aa_ref = shift;
	my %most_freq_ic = ();
	my %most_freq_ic_value = ();

	foreach my $codon (keys %{$icu_ref}) {
		my $aa = $c2aa_ref->{$codon};
		if (! exists $most_freq_ic_value{$aa}) {
			$most_freq_ic_value{$aa} = 0;
		}
		if ($most_freq_ic_value{$aa} < $icu_ref->{$codon}) {
			$most_freq_ic{$aa} = $codon;
			$most_freq_ic_value{$aa} = $icu_ref->{$codon};
		}
	}
	return \%most_freq_ic, \%most_freq_ic_value;
}

sub comp_cc_counts
{
	my $seq = shift;
	my $header = shift;
	my $c2aa_ref = shift;
	my $aa2c_ref = shift;
	my %cc_counts = ();
	my %aap_counts = ();

	foreach my $codon1 (keys %{$c2aa_ref}) {
		if ($c2aa_ref->{$codon1} eq "*") {
			next;
		}
		foreach my $codon2 (keys %{$c2aa_ref}) {
			$cc_counts{$codon1.$codon2} = 0;
		}
	}
	foreach my $aa1 (keys %{$aa2c_ref}) {
		if ($aa1 eq "*") {
			next;
		}
		foreach my $aa2 (keys %{$aa2c_ref}) {
			$aap_counts{$aa1.$aa2} = 0;
		}
	}

	my @codons_array = unpack '(a3)*', $seq;
	my $codon = shift @codons_array;
	if (exists $c2aa_ref->{$codon} && exists $c2aa_ref->{$codons_array[0]} && $c2aa_ref->{$codon} eq "M") {
		++$cc_counts{$codon.$codons_array[0]};
		++$aap_counts{$c2aa_ref->{$codon}.$c2aa_ref->{$codons_array[0]}};
	}
	for (my $i=0; $i<$#codons_array; ++$i) {
		if (exists $c2aa_ref->{$codons_array[$i]} && exists $c2aa_ref->{$codons_array[$i+1]}) {
			if ($c2aa_ref->{$codons_array[$i]} eq "*") {
				next;
			}
			if ($c2aa_ref->{$codons_array[$i+1]} eq "*" && $i+1<$#codons_array) {
				next;
			}
			++$cc_counts{$codons_array[$i].$codons_array[$i+1]};
			++$aap_counts{$c2aa_ref->{$codons_array[$i]}.$c2aa_ref->{$codons_array[$i+1]}};
		}
	}
	return \%cc_counts, \%aap_counts;
}

sub comp_icu
{
	my $c2aa_ref = shift;
	my $icu_ref = shift;
	my $num_ic = shift;
	my $c_counts_ref = shift;
	my $aa_counts_ref = shift;

	my $score = 0.0;
	foreach my $codon (keys %{$icu_ref}) {
		if ($aa_counts_ref->{$c2aa_ref->{$codon}} != 0) {
			$score += abs( ($c_counts_ref->{$codon}/$aa_counts_ref->{$c2aa_ref->{$codon}})
					- $icu_ref->{$codon} );
		} else {
			$score += $icu_ref->{$codon};
		}
	}
	return $score / $num_ic;
}

sub comp_cc
{
	my $c2aa_ref = shift;
	my $cc_ref = shift;
	my $num_cc = shift;
	my $cc_counts_ref = shift;
	my $aap_counts_ref = shift;

	my $score = 0.0;
	foreach my $codonp (keys %{$cc_ref}) {
		if ($aap_counts_ref->{ $c2aa_ref->{substr($codonp,0,3)}.$c2aa_ref->{substr($codonp,3)} } != 0) {
			$score += abs( ($cc_counts_ref->{$codonp}/$aap_counts_ref->{$c2aa_ref->{substr($codonp,0,3)}.$c2aa_ref->{substr($codonp,3)}})
					- $cc_ref->{$codonp} );
		} else {
			$score += $cc_ref->{$codonp};
		}
	}
	return $score / $num_cc;
}

sub comp_cai
{
	my $seq = shift;
	my $c2aa_ref = shift;
	my $icu_ref = shift;
	my $most_freq_ic_value_ref = shift;

	my @codons_array = unpack "(a3)*", $seq;
	my $score = 0.0;
	foreach my $codon (@codons_array) {
		$score += log($icu_ref->{$codon} / $most_freq_ic_value_ref->{$c2aa_ref->{$codon}});
	}
	return exp(1.0 / ($#codons_array+1) * $score);
}

sub comp_hidden
{
	my $seq = shift;
	my $aa2c_ref = shift;

	my $count = 0;
	my @stoplist = @{$aa2c_ref->{"*"}};
	foreach my $codon(@stoplist) {
		for (my $i=0; $i<(length($seq)-3); ++$i) {
			if( substr($seq,$i,3) eq $codon ) {
				++$count;
			}
		}
	}
	return $count;
}

sub comp_gc
{
	my $seq = shift;

	my @array = unpack "(a)*", $seq;
	my $GCcount = 0;
	foreach my $nuc (@array) {
		if ($nuc eq "G" || $nuc eq "C") {
			++$GCcount;
		}
	}
	return 100.0 * $GCcount / ($#array+1);
}

sub comp_exc {
	my $seq = shift;
	my $exclu_ref = shift;

	my $score = 0;
	foreach my $exclseq (@{$exclu_ref}) {
	    my $pos = 0;
		my $matches = 0;
		while (1) {
			$pos = index($seq, $exclseq, $pos);
			last if($pos < 0);
			++$score;
			++$pos;
		}
	}
	return $score;
}

sub comp_repeat
{
	my $seq = shift;
	my $repeat_ref = shift;

	my $score = 0;
	foreach my $rep_len (keys %{$repeat_ref}) {
		my $rep_num = $repeat_ref->{$rep_len};
		for (my $i=0; $i<=(length($seq)-$rep_len*$rep_num); ++$i) {
			if (substr($seq,$i,$rep_len*$rep_num) eq 
				substr($seq,$i,$rep_len) x $rep_num) {
				++$score;
			}
		}
	}
	return $score;
}
























1
