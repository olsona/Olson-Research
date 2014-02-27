#!/usr/bin/perl
my $input = $ARGV[0];
my $prename = substr($input, 0, -7);
my $name = $(split(/\\/, $prename))[-1];
printf $name . "/n";