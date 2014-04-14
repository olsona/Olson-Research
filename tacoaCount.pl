#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;

# Disable buffering of STDOUT
$| = 1;

our $KmerSize = 5;

####################################################################
### Init functions

sub parseArgs {
	my $opts = {};
	getopts('k:', $opts);
	$KmerSize = $opts->{'k'} if defined $opts->{'k'};
}

sub usage {
	print <<EOFUSAGE;
Usage: countKmer.pl [-k 5] <indexFile> <resultsFile>
    For each possible words in the kmer of length -k count the number of times they are found in the fasta sequence file
    -k <size>   size of the kmer to analyze. Default 5
    
EOFUSAGE
	exit(1);
}


####################################################################
### Processing functions

sub processFile {
	my $file = shift;
    
	#generate all possible word for a k-mer
    
	my $k0 = $KmerSize;
	
	my $k0Counts = kmer_generator($k0);
    
	#print STDOUT "Counting the number of $KmerSize nt long kmers in $file\n";
	open(my $FH, $file) || die "Can't open file $file\n";
	my $inSeq = 0;
	my $seq = '';
    my $ct = 0;
    my $gc = 0;
	while(<$FH>) {
		# Clean up newlines
		s/[\r\n]+$//;
        
		if (/^>/) {
			# New sequence
            
			# Process previous sequence, if exists
			if (length($seq) > 0) {
                
                $ct += length($seq);
                
                #Alex's fast fix
                $k0Counts->{$_}++ foreach unpack('(A' . $k0 . 'X' . ($k0 - 1) . ')*', $seq);
                
                $gc += ( $seq =~ tr/CG// );
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
        $ct += length($seq);
        
        #Alex's fast fix
        $k0Counts->{$_}++ foreach unpack('(A' . $k0 . 'X' . ($k0 - 1) . ')*', $seq);
        $gc += ( $seq =~ tr/CG// );
	}
    
	close($FH);
    
	# Return results
	return ($k0Counts, $ct, $gc);
}

sub reverse_complement {
    #http://code.izzid.com/2011/08/25/How-to-reverse-complement-a-DNA-sequence-in-perl.html
    my $dna = shift;
    
	# reverse the DNA sequence
    my $revcomp = reverse($dna);
    
	# complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub countKmer {
    my ($kmerSize, $seq, $kmerCounts) = @_;
	my $beenThere = {};
    
	# Iterate through all kmers in the sequence
	for (my $i = 0; $i <= (length($seq) - $kmerSize); $i++) {
		# Extract the next kmer from the sequence
		my $word = substr($seq, $i, $kmerSize);
        my $rc = reverse_complement($word);
        
		# Skip unless this word is one of our enumerated kmers
		next unless exists $kmerCounts->{$word};
        
        # Count all instances of a kmer per sequence
        # Increment this kmer's count
        $kmerCounts->{$word}++;
        $kmerCounts->{$rc}++;
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
            print STDOUT "WARN: Invalid line format in index: $_\n";
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

# Read each input file
foreach my $inputName (sort keys %$inputFiles) {
	my $inputFile = $inputFiles->{$inputName};
	print STDOUT "Processing File: $inputName - $inputFile\n";
    
	# Read the input file
	my ($kCounts, $s, $gc) = processFile($inputFile);
    
    my $gcpr = $gc/(2*$s);
    
    my %probs = ('C', $gcpr, 'G', $gcpr, 'A', 1-$gcpr, 'T', 1-$gcpr);
    
	# Write results
	print $resultsfh $inputName;
	my @results = ();
	foreach my $o (@allKmers) {
        my $go = 0;
        my $Oo = $kCounts->{$o};
        if ($Oo > 0) {
            my @strarray = unpack 'C*', $o;
            my $pr = 1.0;
            foreach my $c (@strarray) {
                $pr *= %probs{$c};
            }
            my $Eo = $pr*$s;
            if ($Oo > $Eo) {
                $go = $Oo/$Eo;
            }
            else {
                $go = $Eo/$Oo;
                $go *= -1;
            }
        }
        printf $resultsfh ":%.8f", $go;
	}
	print $resultsfh "\n";
}

close($resultsfh);