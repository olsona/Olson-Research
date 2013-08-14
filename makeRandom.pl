#!/usr/bin/perl
my ($p, $input, $db, $seq) = @ARGV;
open(IN, $input);
open(OUT1, '>', $db);
open(OUT2, '>', $seq);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    if ( rand(1) < $p ) {
        print OUT1 ">$name\n";
        print OUT1 "$seq\n";
    } else {
        print OUT2 ">$name\n";
        print OUT2 "$seq\n";
    }
}
close(IN);
close(OUT1);
close(OUT2);