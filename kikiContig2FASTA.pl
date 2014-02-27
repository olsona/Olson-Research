#!/usr/bin/perl
my $input = $ARGV[0];
my $prename = substr($input, 0, -7);
open(IN, $input);
while (<IN>) {
    chomp;
    ($a,$b,$c,$d,$e,$f,$g) = split("\t");
    printf $f . "\n";
}