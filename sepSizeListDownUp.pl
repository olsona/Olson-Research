#!/usr/bin/perl
my ($threshold, $path, $input, $smaller, $bigger) = @ARGV;
open(IN, $input);
open(SMLR1, '>', "$smaller-1");
open(SMLR2, '>', "$smaller-2");
open(BGGR, '>', $bigger);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    $len = length($seq);
    if ( $len >= $threshold ) {
        print BGGR "$name\t$seq\n";
    } else {
        $nname = substr $name, 1;
        $allname = $nname."$len";
        open(FI, '>', "$path$allname.fna");
        print FI "$name\n$seq";
        close(FI);
        print SMLR1 "$path$allname.fna\n";
        print SMLR2 "$allname\t$path$allname.fna\n";
    }
}
close(IN);
close(SMLR1);
close(SMLR2);
close(BGGR);
