#!/usr/bin/perl
my ($input, $output) = @ARGV;
printf $output . "\n";
my $name = substr($input, 0, -7);
$ct = 0;
open(IN, $input);
open(OUT, $output);
while (<IN>) {
    chomp;
    ($a,$sz,$c,$d,$e,$f,$g) = split("\t");
    print OUT "$name_$sz_$ct\n$g\n";
    $ct++;
}
close(IN);
close(OUT);