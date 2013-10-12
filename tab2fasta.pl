#!/usr/bin/perl
my ($input, $output) = @ARGV;
open(IN, $input);
open(OUT, '>', $output);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    print OUT "$name\n";
    print OUT "$seq\n";
}