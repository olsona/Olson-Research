#!/usr/bin/perl

use strict;
use warnings;

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

# Reads counts from jellyfish count files
sub getCounts {
    my ($inFile, $kmerSize) = @_;
    my $kmrHash = kmer_generator($kmerSize);
    open(IN, $inFile);
    while(<IN>) {
        my ($kmr, $ct) = split(" ");
        $kmrHash->{$kmr} = int($ct);
    }
    return $kmrHash;
}


#---MAIN CODE---#
my ($jellyfishPath, $indexFile, $resultsFile) = @ARGV;

# Read index file
my $inputFiles = readIndex($indexFile);

# Open results file
open(my $resultsfh, '>', $resultsFile) or die "FATAL: Unable to write results file: $!\n";

# Main loop
foreach my $inputName (sort keys %$inputFiles) {
    my $inputFile = $inputFiles->{$inputName};
    print STDOUT "Processing File: $inputName - $inputFile\n";
    
    # prepare for running Jellyfish
    my $filesize = -s $inputFile;
    my $hsh = int($filesize/8);
    my $jf4 = "${inputFile}_JF_4";
    my $jf3 = "${inputFile}_JF_3";
    my $jf2 = "${inputFile}_JF_2";
    my $ct4 = "${inputFile}_ct_4";
    my $ct3 = "${inputFile}_ct_3";
    my $ct2 = "${inputFile}_ct_2";
    
    # Run Jellyfish
    system("${jellyfishPath}/jellyfish count -m 4 -s ${hsh} -t 10 -o ${jf4} -C ${inputFile}");
    system("${jellyfishPath}/jellyfish dump -c ${jf4} > ${ct4}");
    system("rm ${jf4}");
    system("${jellyfishPath}/jellyfish count -m 3 -s ${hsh} -t 10 -o ${jf3} -C ${inputFile}");
    system("${jellyfishPath}/jellyfish dump -c ${jf3} > ${ct3}");
    system("rm ${jf3}");
    system("${jellyfishPath}/jellyfish count -m 2 -s ${hsh} -t 10 -o ${jf2} -C ${inputFile}");
    system("${jellyfishPath}/jellyfish dump -c ${jf2} > ${ct2}");
    system("rm ${jf2}");
    
    # Get results from Jellyfish
    my $kmr4 = getCounts($ct4,4);
    my $kmr3 = getCounts($ct3,3);
    my $kmr2 = getCounts($ct2,2);
    
    #for (sort(keys %$kmr4)) {
    #    print $_ . "," . $kmr4->{$_} . "\n";
    #}
    
    my @all4mers = sort keys kmer_generator(4);
    
    # Generate and write results
    print $resultsfh $inputName;
    foreach my $w1h (@all4mers) {
        my $w1h1	= substr($w1h, 0, 3);
        my $w2h1	= substr($w1h, 1, 2);
		my $w2h		= substr($w1h, 1, 3);
        
        my $Nw		= $kmr4->{$w1h};
		my $Nw1h1	= $kmr3->{$w1h1};
		my $Nw2h	= $kmr3->{$w2h};
		my $Nw2h1	= $kmr2->{$w2h1};
        
        my $Nhathat = 0;
        my $Vhathat = 0;
        if ($Nw2h1 >= 1) {
            $Nhathat = ($Nw1h1*$Nw2h)/$Nw2h1;
            $Vhathat = (($Nw1h1*$Nw2h)/($Nw2h1*$Nw2h1*$Nw2h1)) * ($Nw2h1-$Nw1h1) * ($Nw2h1-$Nw2h);
        }
        
        my $ZM      = 0;
        if ($Nhathat > 0 && $Vhathat > 0) {
            $ZM     = ($Nw-$Nhathat)/sqrt($Vhathat);
        }
        
        printf $resultsfh ":%.8f", $ZM;
    }
    print $resultsfh "\n";
}

close($resultsfh);