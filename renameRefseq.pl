#!/usr/bin/perl
my ($inpath, $outpath, $out) = @ARGV;
my @files = <$inpath/*.fna>;
open($OUT1, '>', "$out-1");
open($OUT2, '>', "$out-2");
for my $input (@files) {
    open($IN, $input);
    my $outfh;
    while (<$IN>) {
        chomp;
        my $cr = substr $_, 0, 1;
        #print $cr . "\n";
        if ($cr eq ">" ) {
            print "hello\n";
            close($outfh) if defined $outfh;
            my $name = substr $_, 1;
            print $name . "\n";
            my $spname = (split '|', $name)[4];
            my $cname = (split ',', $spname)[0];
            open($outfh, '>', "$outpath$cname.fna");
            print $outfh ">$cname\t";
            print $OUT1 "$outpath$cname.fna\n";
            print $OUT2 "$cname\t$outpath$cname.fna\n";
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