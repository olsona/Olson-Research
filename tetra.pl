#!/usr/bin/perl

use strict;
use warnings;
3use Statistics::Basic qw(:all nofill);

my ($DBfile) = @ARGV;

print $DBfile;

my @DBlines = `cat $DBfile`;
my @DBvectors = ();
foreach my $dbl (@DBlines) {
    #print $dbl;
    my @dbln = split ':', $dbl;
    shift @dbln;
    my $dbv = vector(@dbln);
    print "$dbv\n";
}