#!/usr/bin/perl
my ($threshold, $path, $input, $smaller, $bigger) = @ARGV;
use Text::Wrap;
open(IN, $input);
open(SMLR, '>', $smaller);
open(BGGR, '>', $bigger);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    if ( length($seq) >= $threshold ) {
        print BGGR "$name\t$seq\n";
    } else {
        $nname = substr $name, 1;
        open(FI, '>', "$path$nname.fna");
        print FI "$name\n";
        $Text::Wrap::columns = 60;
        print FI wrap('','',$seq);
        close(FI);
        print SMLR "$nname\t$path$nname.fna\n";
    }
}
close(IN);
close(SMLR);
close(BGGR);