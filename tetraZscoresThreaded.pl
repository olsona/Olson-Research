#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;
use threads::shared;

# Disable buffering of STDOUT
$| = 1;

our $NumJobs = 7;

our $KmerSize = 6;
our $OnlyCountFirstKmerOccurrence = 1;



####################################################################
### Init functions

sub parseArgs {
	my $opts = {};
	getopts('j:k:m', $opts);
	$NumJobs = $opts->{'j'} if defined $opts->{'j'};
	$KmerSize = $opts->{'k'} if defined $opts->{'k'};
	$OnlyCountFirstKmerOccurrence = 0 if defined $opts->{'m'};
}

sub usage {
	print <<EOFUSAGE;
Usage: countKmer.pl [-k 6 -m] <indexFile> <resultsFile>
For each possible words in the kmer of length -k count the number of time they are found in the fasta sequence file
  -j <jobs>   number of threads to use
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
	#printThread("Generating all the posible sequences $KmerSize long\n");
	
	my $k0 = $KmerSize;
	my $k1 = $KmerSize - 1;
	my $k2 = $KmerSize - 2;
	
	my $k0Counts = kmer_generator($k0);
	my $k1Counts = kmer_generator($k1);
	my $k2Counts = kmer_generator($k2);

	#printThread("Counting the number of $KmerSize nt long kmers in $file\n");
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
	my $max = length($seq) - $kmerSize;
	my $word;

	if ($OnlyCountFirstKmerOccurrence) {
		# Count only one occurrence of a kmer per sequence
		my $beenThere = {};

		# Iterate through all kmers in the sequence
		for (my $i = 0; $i <= $max; $i++) {
			# Extract the next kmer from the sequence
			$word = substr($seq, $i, $kmerSize);

			if ( !exists $beenThere->{$word} ) {
				# Increment this kmer's count
				$kmerCounts->{$word}++;
				# Mark this kmer as seen
				$beenThere->{$word} = 1;
			}
		}

	} else {
		# Count all instances of a kmer per sequence

		# Iterate through all kmers in the sequence
		for (my $i = 0; $i <= $max; $i++) {
			# Extract the next kmer from the sequence
			$word = substr($seq, $i, $kmerSize);

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

my $PrintLock :shared;
sub printThread {
	lock($PrintLock);
	print STDERR "tid=" . threads->tid() . " ", @_;
}



####################################################################
### Main code

# Read command-line args
parseArgs();
my ($indexFile, $resultsFile) = @ARGV;
die usage() unless -f $indexFile && $resultsFile ne '';

# Read index file
my $inputFiles = readIndex($indexFile);

# Generate list of all kmers for results printing
our @AllK0Kmers = sort keys kmer_generator($KmerSize);
our @AllK1Kmers = sort keys kmer_generator($KmerSize - 1);
our @AllK2Kmers = sort keys kmer_generator($KmerSize - 2);

our $h_1 = $KmerSize-1;

my $inQueue = Thread::Queue->new();
my $outQueue = Thread::Queue->new();

# Spawn worker threads
my $threads = {};
foreach my $tid (1..$NumJobs) {
	
	$threads->{$tid} = threads->create( sub {
		# Continually get new data to process from inQueue
		while ( my $work = $inQueue->dequeue() ) {
			my ($inputName, $inputFile) = @$work;
			printThread("Processing File: $inputName - $inputFile\n");

			# Read the input file
			my ($k0Counts, $k1Counts, $k2Counts) = processFile($inputFile);
			
			# Write results
			my $results = $inputName;
			foreach my $kmer (@AllK0Kmers) {
				# I have no idea what's going on here, what the intention is, etc...
				my $w1h_1	= substr($kmer, 0, $h_1);
				my $w2h_1	= substr($kmer, 1, $h_1-1);
				my $w2h		= substr($kmer, 1, $h_1);
				my $Nw		= $k0Counts->{$kmer};
				my $Nw1h1	= $k1Counts->{$w1h_1};
				my $Nw2h	= $k1Counts->{$w2h};
				my $Nw2h1	= $k2Counts->{$w2h_1};
				#printf "%d %d %d %d\n", $Nw, $Nw1h1, $Nw2h, $Nw2h1;
				$results .= sprintf(":%.d", $Nw);
			}

			# Store result
			$outQueue->enqueue($results);
		}
	});
}


# Read each input file
foreach my $inputName (sort keys %$inputFiles) {
	my $inputFile = $inputFiles->{$inputName};

	# Enqueue each input file
	$inQueue->enqueue( [$inputName, $inputFile] );
}
# Indicate that the inQueue is complete
$inQueue->end();

# Open results file
open(my $resultsfh, '>', $resultsFile) or die "FATAL: Unable to write results file: $!\n";

# Check for results until all threads are done and all queues are empty
while (1) {
	# Check if all threads are done
	my $runningThreads = threads->list(threads::running);
	last if $runningThreads <= 0 && (!defined $inQueue->pending() || $inQueue->pending() <= 0) && $outQueue->pending() <= 0;

	# Poll for a result every 1sec
	my $result = $outQueue->dequeue_timed(1);
	if (defined $result) {
		# Print result to file
		print $resultsfh $result, "\n";
	}
}

# Join all threads
foreach my $thread ( threads->list() ) {
	$thread->join();
}

close($resultsfh);


