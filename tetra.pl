#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Basic qw(:all nofill);

my ($DBfile) = @ARGV;

my @DBlines = `cat $DBfile`;
my @DBvectors = ();
my @DBnames = ();
foreach my $dbl (@DBlines) {
    #print $dbl;
    my @dbln = split ':', $dbl;
    push(@DBnames, $dbln[0]);
    shift @dbln;
    my $dbv = vector(@dbln);
    push(@DBvectors, $dbv);
}

print "@DBnames\n";