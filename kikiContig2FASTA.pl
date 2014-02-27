#!/usr/bin/perl
my $input = $ARGV[0];
my $prename = substr($input, 0, -7);
printf $prename;