#!/usr/bin/perl
my ($threshold, $path, $input, $smaller, $bigger) = @ARGV;
open(IN, $input);
open(BGGR1, '>', "$bigger-1");
open(BGGR2, '>', "$bigger-2");
open(SMLR, '>', $smaller);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    if ( length($seq) >= $threshold ) {
        print SMLR "$name\t$seq\n";
    } else {
        $nname = substr $name, 1;
        open(FI, '>', "$path$nname.fna");
        print FI "$name\n$seq";
        close(FI);
        print BGGR1 "$path$nname.fna\n";
        print BGGR2 "$nname\t$path$nname.fna\n";
    }
}
close(IN);
close(BGGR1);
close(BGGR2);
close(SMLR);
