#!/usr/bin/perl

use strict;
use warnings;

my ($infile, $k, $tcount, $countfile) = @ARGV;

open (IN, $infile);

my $c = 0;
my $g = 0;
my $total = 0;
while (<IN>) {
    chomp;
    my $_ = lc($_);
    $total += length;
    $c += ($_ =~ tr/c//);
    $g += ($_ =~ tr/g//);
}

close(IN);
my $gc = $g+$c;

