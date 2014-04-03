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


#---MAIN CODE---#
my ($jellyfishPath, $indexFile, $workingPrefix, $resultsFile) = @ARGV;

my $inputFiles = readIndex($indexFile);

foreach my $inputName (sort keys %$inputFiles) {
    my $inputFile = $inputFiles->{$inputName};
    my $filesize = -s $inputFile;
    my $hsh4 = int($filesize/4);
    my $work4 = "${workingPrefix}_${inputFile}_4";
    echo $work4;
    system("$jellyfishPath/jellyfish count -m 4 -s $hsh4 -t 10 -o $work4 -C $inputFile");
}