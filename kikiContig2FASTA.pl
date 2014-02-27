#!/usr/bin/perl
my ($input, $output) = @ARGV;
open(IN, $input);
open(OUT, '>', $output);
my $name = substr($input, 0, -7);
my $ct = 0;
while (<IN>) {
    chomp;
    ($a,$sz,$c,$d,$e,$f,$g) = split("\t");
    print OUT "$name_$sz_$ct\n$g\n";
    $ct++;
}
close(IN);
close(OUT);