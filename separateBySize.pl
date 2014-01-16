#!/usr/bin/perl
my ($input, $input, $smaller, $bigger) = @ARGV;
open(IN, $input);
open(OUT1, '>', $smaller);
open(OUT2, '>', $bigger);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    if ( length($seq) >= $threshold ) {
        print OUT2 "$name\t$seq\n";
    } else {
        print OUT1 "$name\t$seq\n";
    }
}
close(IN);
close(OUT1);
close(OUT2);