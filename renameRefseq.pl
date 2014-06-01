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
            my $spname = (split $d, $name)[4];
            my $ge = (split " ", $spname)[0];
            my $sp = (split " ", $spname)[1];
            if ($sp eq ".sp") {
                $sp = (split " ", $spname[2];
            }
            my $fname = (substr $ge, 1)."_".$sp;
            print $fname . "\n";
            if (-e $outfh) {
                open($outfh, '>>', "$outpath$fname.fna");
            }
            else {
                open($outfh, '>', "$outpath$fname.fna");
                print $outfh ">$fname\t";
            }
            print $OUT1 "$outpath$fname.fna\n";
            print $OUT2 "$fname\t$outpath$fname.fna\n";
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