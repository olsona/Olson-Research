#!/usr/bin/perl

use strict;
use warnings;
use File::Slurp;

my ($DBfile) = @ARGV;

print $DBfile;

my @DBlines = read_file($DBfile, chomp => 1);
my @DBvectors = ();
foreach my $dbl (@DBlines) {
    #print $dbl;
    my @dbln = split ':', $dbl;
    shift @dbln;
    my @dbv = vector(@dbln);
}