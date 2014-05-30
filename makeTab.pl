#!/usr/bin/perl
my ($inpath, $outpath, $out) = @ARGV;
my @files = <$inpath/*.fna>;
open($OUT1, '>', "$out-1");
open($OUT2, '>', "$out-2");
my $d = qr/\|/;

for my $input (@files) {
    open($IN, $input);
    my $outfh;
    while (<$IN>) {
        chomp;
        my $cr = substr $_, 0, 1;
        #print $cr . "\n";
        if ($cr eq ">" ) {
            close($outfh) if defined $outfh;
            my $name = substr $_, 1;
            print $name . "\n";
            open($outfh, '>', "$outpath$name.fna");
            print $outfh ">$name\t";
            print $OUT1 "$outpath$name.fna\n";
            print $OUT2 "$name\t$outpath$name.fna\n";
        }
        else {
            #print "hi\n";
            print $outfh $_;
        }
    }
}
close($outfh) if defined $outfh;
close($IN);
close($OUT1);
close($OUT2);