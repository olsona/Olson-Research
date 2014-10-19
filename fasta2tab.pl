#!/usr/bin/perl
my ($input, $output) = @ARGV;
open(IN, $input);
open(OUT, '>', $output);
my $newOrg = 0;
my $name = "";
my $seq = "";
while (<IN>) {
    chomp;
    if ($newOrg == 0) {
        $name = $_;
        $newOrg = 1;
    } else {
        $seq = $_;
        $newOrg = 0;
        print OUT "$name\t$seq\n";
    }
}
close(IN);
close(OUT);