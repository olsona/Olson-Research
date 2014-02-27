#!/usr/bin/perl
my $input = @ARGV;
my $name = substr($input, 0, -3);
printf $name;