#!/usr/bin/perl
my $input = $ARGV[0];
my $name = split(substr($input, 0, -7),"/")[-1];
printf $name . "/n";