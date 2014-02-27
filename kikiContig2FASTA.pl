#!/usr/bin/perl
my $input = @ARGV;
open(IN, $input);
name=substr($input, 0, -7);
printf $name . "\n";