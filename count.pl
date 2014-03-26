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
    
    my $k = $KmerSize;
    my $k1 = $KmerSize-1;
    my $k2 = $KmerSize-2;
    
	my $kmerCounts = kmer_generator($k);
    my $kminus1Counts = kmer_generator($k1);
    my $kminus2Counts = kmer_generator($k2);

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
				countKmer(\$seq, $kmerCounts);
                countKmer(\$seq, $kminus1Counts);
                countKmer(\$seq, $kminus2Counts);
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
		countKmer(\$seq, $kmerCounts);
        countKmer(\$seq, $kminus1Counts);
        countKmer(\$seq, $kminus2Counts);
	}

	close($FH);

    my %allKmerInfo = ($k, $kmerCounts, $k1, $kminus1Counts, $k2, $kminus2Counts);
    
	# Return results
	return (%allKmerInfo);
}

sub countKmer {
	my $seq = ${ $_[0] };
	my $kmerCounts = $_[1];
	my $kmerLength = $KmerSize;
	my $beenThere = {};

	# Iterate through all kmers in the sequence
	for (my $i = 0; $i <= (length($seq) - $kmerLength); $i++) {
		# Extract the next kmer from the sequence
		my $word = substr($seq, $i, $kmerLength);

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
my @allKmers = sort keys kmer_generator($KmerSize);
my @minus1Kmers = sort keys kmer_generator($KmerSize - 1);
my @minus2Kmers = sort keys kmer_generator($KmerSize - 2);

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
	my (%kmerInfo) = processFile($inputFile);

    my $kmerCounts = $kmerInfo{$KmerSize};
    my $kminus1Counts = $kmerInfo{$KmerSize-1};
    my $kminus2Counts = $kmerInfo{$KmerSize-2};
    
	# Write results
	print $resultsfh $inputName;
	my @results = ();
	foreach my $kmer (@allKmers) {
        my $w1h_1 = substr $kmer, 0, $h_1;
        my $w2h_1 = substr $kmer, 1, $h_1-1;
        my $w2h = substr $kmer, 1, $h_1;
		my $Nw = $kmerCounts->{$kmer};
        my $Nw1h1 = $kminus1Counts->{$w1h_1};
        my $Nw2h = $kminus1Counts->{$w2h};
        my $Nw2h1 = $kminus2Counts->{$w2h_1};
        printf "%d %d %d %d\n", $Nw, $Nw1h1, $Nw2h, $Nw2h1;
		printf $resultsfh ":%.d", $Nw;
	}
	print $resultsfh "\n";

}

close($resultsfh);


