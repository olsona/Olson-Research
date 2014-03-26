#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Basic qw(:all nofill);

my ($DBfile,$matchfile,$outfile) = @ARGV;

open (OUT, '>', $outfile);

my @DBlines = `cat $DBfile`;
my @DBvectors = ();
my @DBnames = ();
foreach my $dbl (@DBlines) {
    my @dbln = split ':', $dbl;
    push(@DBnames, $dbln[0]);
    shift @dbln;
    my $dbv = vector(@dbln);
    push(@DBvectors, $dbv);
}

my @matchlines = `cat $matchfile`;
my @matchnames = ();
my @matchvectors = ();
foreach my $mal (@matchlines) {
    my @maln = split ':', $mal;
    push(@matchnames, $maln[0]);
    shift @maln;
    my $mav = vector(@maln);
    push(@matchvectors, $mav);
}

my $dbstr = join(',',@DBnames);
print OUT "$dbstr\n";
my $mstr = join(',',@matchnames);
print OUT "$mstr\n";

my $dbLen = scalar $dbVectors;

foreach my $mv (@matchvectors) {
    my %zrec = ();
    for (my $i=0; $i < dbLen; i++) {
        $dv = $dbVectors[i];
        my $cor = corr($mv, $dv);
        my $zc = $cor->query;
        $zrec->{$zc} = $i;
    }
    my @zsort = (sort {$zrec{$a} <=> $zrec{$b}} keys %zrec);
    my $zstr = join(',',@zres);
    print OUT "$zstr\n";
}

close(OUT);