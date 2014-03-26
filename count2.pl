#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;

# Disable buffering of STDOUT
$| = 1;

our $KmerSize = 6;
our $OnlyCountFirstKmerOccurrence = 1;



####################################################################
### Init functions

sub parseArgs {
	my $opts = {};
	getopts('k:m', $opts);
	$KmerSize = $opts->{'k'} if defined $opts->{'k'};
	$OnlyCountFirstKmerOccurrence = 0 if defined $opts->{'m'};
}

sub usage {
	print <<EOFUSAGE;
Usage: countKmer.pl [-k 6 -m] <indexFile> <resultsFile>
For each possible words in the kmer of length -k count the number of time they are found in the fasta sequence file
  -k <size>   size of the kmer to analyze. Default 6
  -m          will count all possible kmer per sequences. Default: only one kmer is counted per sequence entries

EOFUSAGE
	exit(1);
}



####################################################################
### Processing functions

sub processFile {
	my $file = shift;

	#generate all possible word for a k-mer
	print STDERR "Generating all the posible sequences $KmerSize long\n";
	
	my $k0 = $KmerSize;
	my $k1 = $KmerSize - 1;
	my $k2 = $KmerSize - 2;
	
	my $k0Counts = kmer_generator($k0);
	my $k1Counts = kmer_generator($k1);
	my $k2Counts = kmer_generator($k2);

	print STDERR "Counting the number of $KmerSize nt long kmers in $file\n";
	open(my $FH, $file) || die "Can't open file $file\n";
	my $inSeq = 0;
	my $seq = '';
	while(<$FH>) {
		# Clean up newlines
		s/[\r\n]+$//;

		if (/^>/) {
			# New sequence

			# Process previous sequence, if exists
			if (length($seq) > 0) {
				countKmer($k0, $seq, $k0Counts);
				countKmer($k1, $seq, $k1Counts);
				countKmer($k2, $seq, $k2Counts);
			}

			# Setup for reading new sequence
			$seq = '';
			$inSeq = 1;

		} elsif ($inSeq) {
			# Append line to sequence buffer
			$seq .= $_;
		}
	}

	# Process final sequence, if exists
	if (length($seq) > 0) {
        countKmer($k0, $seq, $k0Counts);
        countKmer($k1, $seq, $k1Counts);
        countKmer($k2, $seq, $k2Counts);
	}

	close($FH);

	# Return results
	return ($k0Counts, $k1Counts, $k2Counts);
}

sub countKmer {
    my ($kmerSize, $seq, $kmerCounts) = @_;
	my $beenThere = {};

	# Iterate through all kmers in the sequence
	for (my $i = 0; $i <= (length($seq) - $kmerSize); $i++) {
		# Extract the next kmer from the sequence
		my $word = substr($seq, $i, $kmerSize);

		# Skip unless this word is one of our enumerated kmers
		next unless exists $kmerCounts->{$word};

		if ($OnlyCountFirstKmerOccurrence) {
			# Count only one occurrence of a kmer per sequence

			if ( !exists $beenThere->{$word} ) {
				# Increment this kmer's count
				$kmerCounts->{$word}++;
				# Mark this kmer as seen
				$beenThere->{$word} = 1;
			}

		} else {
			# Count all instances of a kmer per sequence

			# Increment this kmer's count
			$kmerCounts->{$word}++;
		}
	}
}

# Generate a hashref of each kmer possibility
sub kmer_generator {
	my $k = shift;
	my $kmer_p;

	my @bases = ('A','C','G','T');
	my @words = @bases;

	for (my $i=1;$i<$k;$i++) {
		my @newwords;
		foreach my $w (@words) {
			foreach my $b (@bases) {
				push (@newwords,$w.$b);
			}
		}
		undef @words;
		@words = @newwords;
	}
	map{ $kmer_p->{$_}=0} @words;
	return($kmer_p);
}

# Read given index file and return a hashref of names and filenames
sub readIndex {
	my ($indexFile) = @_;
	my $inputFiles = {};
	open(my $fh, $indexFile) or die;
		while (<$fh>) {
			s/[\r\n]+//;
			# Line Format: <Name>\t<Filename>
			if (/^(.+?)\t+(.+)$/) {
				my ($inputName, $inputFile) = ($1, $2);
				die "FATAL: Files with duplicate names in index\n" if exists $inputFiles->{$inputName};
				$inputFiles->{$inputName} = $inputFile;
			} else {
				print STDERR "WARN: Invalid line format in index: $_\n";
			}
		}
	close($fh);
	return $inputFiles;
}



####################################################################
### Main code

# Read command-line args
parseArgs();
my ($indexFile, $resultsFile) = @ARGV;
die usage() unless -f $indexFile && $resultsFile ne '';

# Read index file
my $inputFiles = readIndex($indexFile);

# Open results file
open(my $resultsfh, '>', $resultsFile) or die "FATAL: Unable to write results file: $!\n";

# Generate list of all kmers for results printing
my @allK0Kmers = sort keys kmer_generator($KmerSize);
my @allK1Kmers = sort keys kmer_generator($KmerSize - 1);
my @allK2Kmers = sort keys kmer_generator($KmerSize - 2);

# Write results header
#print $resultsfh 'Name';
#my @results = ();
#foreach my $kmer (@allKmers) {
#	print $resultsfh ':' . $kmer;
#}
#print $resultsfh "\n";

my $h_1 = $KmerSize-1;

# Read each input file
foreach my $inputName (sort keys %$inputFiles) {
	my $inputFile = $inputFiles->{$inputName};
	print STDERR "Processing File: $inputName - $inputFile\n";

	# Read the input file
	my ($k0Counts, $k1Counts, $k2Counts) = processFile($inputFile);
	
	# Write results
	print $resultsfh $inputName;
	my @results = ();
	foreach my $w1h (@allK0Kmers) {
	    # I have no idea what's going on here, what the intention is, etc...
		my $w1h1	= substr($w1h, 0, $h_1);
		my $w2h1	= substr($w1h, 1, $h_1-1);
		my $w2h		= substr($w1h, 1, $h_1);
        
		my $Nw		= $k0Counts->{$w1h};
		my $Nw1h1	= $k1Counts->{$w1h1};
		my $Nw2h	= $k1Counts->{$w2h};
		my $Nw2h1	= $k2Counts->{$w2h1};
    
        my $Nhathat = ($Nw1h1*$Nw2h)/$Nw2h1;
        
		printf $resultsfh ":%.10f", $Nhathat;
	}
	print $resultsfh "\n";

}

close($resultsfh);

