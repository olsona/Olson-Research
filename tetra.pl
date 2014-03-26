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

my $dbLen = @DBvectors;

foreach my $mv (@matchvectors) {
    my %zrec = ();
    foreach my $i (0..$dbLen-1) {
        my $dv = $DBvectors[$i];
        my $cor = corr($mv, $dv);
        my $zc = $cor->query;
        $zrec{$zc} = $i;
    }
    my @zsort = (sort {$b <=> $a} keys %zrec);
    my @zfull = ();
    foreach my $z (@zsort) {
        my $str = sprintf("%.8f:%d",$z,$zrec{$z});
        push(@zfull,$str);
    }
    my $zstr = join(', ',@zfull);
    print OUT "$zstr\n";
}


close(OUT);