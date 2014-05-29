#!/usr/bin/perl
my ($inpath, $outpath, $out) = @ARGV;
my @files = <$inpath/*.fna>;
open($OUT1, '>', "$out-1");
open($OUT2, '>', "$out-2");
for my $input (@files) {
    open($IN, $input);
    my $FI;
    while (<$IN>) {
        chomp;
        my $chr = substr $_, 0, 1;
        if ( $chr == '>' ) {
            close($FI) if (defined $FI);
            my $name = substr $_, 1;
            my $spname = (split '|', $name)[4];
            my $cname = (split ',', $spname)[0];
            open($FI, '>', "$outpath$cname.fna");
            print $FI ">$cname\t";
            print $OUT1 "$outpath$cname.fna\n";
            print $OUT2 "$cname\t$outpath$cname.fna\n";
        }
        else {
            print $FI "$_";
        }
    }
}
close($IN);
close($OUT1);
close($OUT2);