#!/usr/bin/perl
my ($path, $input, $out) = @ARGV;
open(IN, $input);
open(OUT1, '>', "$out-1");
open(OUT2, '>', "$out-2");
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    $nname = substr $name, 1;
    $spname = split('_len',$nname)[0];
    open(FI, '>', "$path$nname.fna");
    print FI "$name\n$seq";
    print $spname."\n";
    close(FI);
    print OUT1 "$path$nname.fna\n";
    print OUT2 "$nname\t$path$nname.fna\n";
}
close(IN);
close(OUT1);
close(OUT2);