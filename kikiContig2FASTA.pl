#!/usr/bin/perl
my $input = $ARGV[0];
my $name = substr($input, 0, -7);
printf $name;