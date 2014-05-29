#!/usr/bin/perl
use File::Slurp;
my ($inpath, $outpath, $out) = @ARGV;
my @files = read_dir $inpath;
open(OUT1, '>', "$out-1");
open(OUT2, '>', "$out-2");
for my $input (@files) {
    open(IN, $input);
    while (<IN>) {
        chomp;
        ($name, $seq) = split("\t");
        $spname = (split '|', $name)[4];
        $cname = (split ',' $spname)[0];
        open(FI, '>', "$outpath$cname.fna");
        print FI ">$cname\n$seq\n";
        close(FI);
        print OUT1 "$outpath$cname.fna\n";
        print OUT2 "$cname\t$outpath$cname.fna\n";
    }
}
close(IN);
close(OUT1);
close(OUT2);