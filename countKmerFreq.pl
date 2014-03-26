#!/usr/bin/perl
# countKmer.pl take a series of fasta files and count all the occurrence of all the
# possible kmers in each files. It's a multi-threaded scripts that treat multiple files 
# asynchronously

use strict;
use warnings;
use Getopt::Std;
use threads;
use Thread::Queue;

our ($opt_k,$opt_m,$opt_p, $opt_f);
our $SEP="\t";

MAIN:{
	init();

	my $q = Thread::Queue->new;
	$q->enqueue(@ARGV);

	my $num_workers = @ARGV < $opt_p ? @ARGV : $opt_p; # no need to be wasteful :)

	for (1 .. $num_workers) {
		threads->new(\&worker, $q);
	}

	$_->join for threads->list;
}

sub worker {
	my $queue = shift;
	while (my $filename = $queue->dequeue_nb) {
		processFile($filename);
	}
	return(0);
}

sub processFile {
	my $file = shift;

	#generate all possible word for a k-mer (define at run time by $opt_k)
	print STDERR "Generating all the posible sequences $opt_k long\n";
	my $kmer_p = kmer_generator($opt_k);

	my $kmerCount = 0;

	print STDERR "Counting the number of $opt_k nt long kmers in $file\n";
	loadSeq($file,$kmer_p, \$kmerCount);

	##print out the hits
	if ($opt_f) {
		print STDERR "Outputting frequencies\n";
	} else {
		print STDERR "Outputting raw counts\n";
	}
	printHits($kmer_p, \$kmerCount, $file);
}

sub loadSeq {
	my $file = shift;
	my $kmer_p = shift;
	my $kmerCount_p = shift;

	open FH, "$file" || die "Can't open file $file\n";
	my $f=0;
	my $seq;
	while(<FH>) {
		if(/^>/) {
			countKmer(\$seq,$kmer_p, $kmerCount_p) if $seq;
			$seq='';
			$f=1;
		}elsif($f==1) {
			chomp;
			next if "";
			$seq = join("",$seq,$_);
		}
	}
	countKmer(\$seq,$kmer_p, $kmerCount_p) if $seq; # do not forget the last sequence
	close FH;
	return(0);
}

sub countKmer {
	my $seq_p = shift;
	my $kmer_p = shift;
	my $kmerCount_p = shift;
	my $k = $opt_k;
	my %beenThere;

	for (my $i=0;$i <= length(${$seq_p})-$k;$i++) {
		my $w = substr(${$seq_p},$i,$k);
		unless ($opt_m) {
			#Count only one occurrence of a kmer per sequence
			if ( !exists $beenThere{$w} && exists $kmer_p->{$w} ) {
				# Increment this kmer's count
				$kmer_p->{$w}++;
				# Increment the total count
				${$kmerCount_p}++;
			}
			$beenThere{$w}=1;
		}else{
			#Count all instances of a kmer per sequence
			if ( exists $kmer_p->{$w} ) {
				# Increment this kmer's count
				$kmer_p->{$w}++;
				# Increment the total count
				${$kmerCount_p}++;
			}
		}
	}
	return(0);
}

sub printHits {
	my $kmer_p = shift;
	my $kmerCount_p = shift;
	my $file = shift;

	##print out the hits
	my ($dir,$pre,$suf) = ($file =~ /(^.+\/|^)(.+)\.(.+$)/);
	#open OUT, ">$pre.hits" || die "Can't create file $pre.hits\n";

	for my $key ( sort keys %{$kmer_p} ) {
		if ($opt_f) {
			# Output frequencies
			my $freq = ( $kmer_p->{$key} / $$kmerCount_p );
			#print OUT join($SEP, $key, $freq), "\n";
            printf "%.10f:", $freq;
		} else {
			# Output raw counts
			#print OUT join($SEP, $key, $kmer_p->{$key}), "\n";
            print $kmer_p->{$key}, ":";
		}
	}

	#close OUT;
    print "\n";

	return(0);
}


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

sub init {
	getopts("p:k:mf");
	unless (@ARGV) {
		print("\nUsage: countKmer.pl [-k 6 -p 4 -m] sequence_1.fa [sequence_2.fa ...]\n",
		"\tFor each possible words in the kmer of lenght -k,\n",
		"\tcount the number of time they are found in the fasta sequence file\n",
		"\t-k\tsize of the kmer to analyze. Default 6\n",
		"\t-f\toutput frequencies instead of raw counts.\n",
		"\t-m\twill count all possible kmer per sequences.\n",
		"\t\tDefault: only one kmer is counted per sequence entries\n",
		"\t-p\tThe number of jobs to process simultaneoulsy. Normally, the number of available processors\n",
		"\t\tDefault: 4\n",
		"\n",
		);
		exit(1)
	}
	$opt_k=6 unless $opt_k;
	$opt_p=4 unless $opt_p;
	return(0);
}

