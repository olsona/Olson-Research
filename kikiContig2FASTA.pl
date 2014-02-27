#!/usr/bin/perl
my $input = @ARGV;
printf $input;
my $name = substr($input, 0, -3);
printf $name;