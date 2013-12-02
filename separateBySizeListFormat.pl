#!/usr/bin/perl
my ($threshold, $path, $input, $smaller, $bigger) = @ARGV;
open(IN, $input);
open(SMLR1, '>', "$smaller-1");
open(SMLR2, '>', "$smaller-2");
open(BGGR, '>', $bigger);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    if ( length($seq) >= $threshold ) {
        print BGGR "$name\t$seq\n";
    } else {
        $nname = substr $name, 1;
        open(FI, '>', "$path$nname.fna");
        print FI "$name\n$seq";
        close(FI);
        print SMLR1 "$path$nname.fna\n";
        print SMLR2 "$nname\t$path$nname.fna\n";
    }
}
close(IN);
close(SMLR1);
close(SMLR2);
close(BGGR);
