#!/usr/bin/perl
my ($threshold, $input, $db, $seq, $nums) = @ARGV;
my ($dbCt, $seqCt) = (0,0);
open(IN, $input);
open(OUT1, '>', $db);
open(OUT2, '>', $seq);
open(OUT3, '>', $nums);
while (<IN>) {
    chomp;
    ($name, $seq) = split("\t");
    if ( length($seq) >= $threshold ) {
        print "Hello";
        print OUT1 "$name\n";
        print OUT1 "$seq\n";
        $dbCt ++;
    } else {
        print OUT2 "$name\n";
        print OUT2 "$seq\n";
        $seqCt ++;
    }
}
print OUT3 "$dbCt,$seqCt\n";
close(IN);
close(OUT1);
close(OUT2);
close(OUT3);