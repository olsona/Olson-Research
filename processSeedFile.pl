#!/usr/bin/perl
my ($path, $input, $output) = @ARGV;
open(IN, $input);
open(OUT2, '>', "$output-2");
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    $nname = substr $name, 1;
    open(FI, '>', "$path$nname.fna");
    print FI "$name\n$seq";
    close(FI);
    print OUT2 "$nname\t$path$nname.fna";
}
close(IN);
close(OUT2);