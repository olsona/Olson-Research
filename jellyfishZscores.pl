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

sub getCounts {
    my ($inFile) = @_;
    my $kmrHash;
    open(IN, $inFile);
    while(<IN>) {
        my ($kmr, $ct) = split(" ");
        $kmrHash->{$kmr} = int($ct);
    }
    return $kmrHash;
}


#---MAIN CODE---#
my ($jellyfishPath, $indexFile, $resultsFile) = @ARGV;

my $inputFiles = readIndex($indexFile);

foreach my $inputName (sort keys %$inputFiles) {
    my $inputFile = $inputFiles->{$inputName};
    my $filesize = -s $inputFile;
    my $hsh = int($filesize/4);
    my $jf4 = "${inputFile}_JF_4";
    my $jf3 = "${inputFile}_JF_3";
    my $jf2 = "${inputFile}_JF_2";
    my $ct4 = "${inputFile}_ct_4";
    my $ct3 = "${inputFile}_ct_3";
    my $ct2 = "${inputFile}_ct_2";
    
    system("${jellyfishPath}/jellyfish count -m 4 -s ${hsh} -t 10 -o ${jf4} -C ${inputFile}");
    system("${jellyfishPath}/jellyfish dump -c ${jf4} > ${ct4}");
    system("rm ${jf4}");
    system("${jellyfishPath}/jellyfish count -m 3 -s ${hsh} -t 10 -o ${jf3} -C ${inputFile}");
    system("${jellyfishPath}/jellyfish dump -c ${jf3} > ${ct3}");
    system("rm ${jf3}");
    system("${jellyfishPath}/jellyfish count -m 2 -s ${hsh} -t 10 -o ${jf2} -C ${inputFile}");
    system("${jellyfishPath}/jellyfish dump -c ${jf2} > ${ct2}");
    system("rm ${jf2}");
    
    my $kmr4 = getCounts($ct4);
    my $kmr3 = getCounts($ct3);
    my $kmr2 = getCounts($ct2);
    
    foreach my $w1h (sort (keys %$kmr4)) {
        print $w1h."\n";
    }
}