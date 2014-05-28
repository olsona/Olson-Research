#!/usr/bin/perl
my ($path, $input, $out) = @ARGV;
open(IN, $input);
open(OUT1, '>', "$out-1");
open(OUT2, '>', "$out-2");
$ct=0;
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    $nname = substr $name, 1;
    $spname = (split '_len', $nname)[0];
    open(FI, '>', "$path$spname.fna");
    print FI "$>spname\n$seq";
    close(FI);
    print OUT1 "$path$spname.fna\n";
    print OUT2 "$spname\t$path$spname.fna\n";
    $ct += 1;
}
close(IN);
close(OUT1);
close(OUT2);
print $ct + "\n";