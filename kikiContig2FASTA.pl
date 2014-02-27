#!/usr/bin/perl
my $input = @ARGV;
my $name = substr($input, 0, -7);
printf $name . "\n";