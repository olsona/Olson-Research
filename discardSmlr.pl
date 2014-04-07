#!/usr/bin/perl
my ($threshold, $path, $input, $bigger) = @ARGV;
open(IN, $input);
open(BGGR, '>', $bigger);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    $len = length($seq);
    if ( $len >= $threshold ) {
        print BGGR "$name\t$seq\n";
    }
}
close(IN);
close(BGGR);
