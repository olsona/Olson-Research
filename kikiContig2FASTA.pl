#!/usr/bin/perl
my $input = $ARGV[0];
my $prename = substr($input, 0, -7);
my @fname = split(/\\/, $prename);
printf $fname . "/n";