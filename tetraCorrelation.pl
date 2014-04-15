#!/usr/bin/perl
use strict;
use warnings;
#use Statistics::Basic qw(:all nofill);
use PDL;
use PDL::Stats;
############################################################
### Functions
# Read DB file input
sub readFile {
	my ($file) = @_;
	open(my $fh, $file);
	my $names	= [];
	my $vectors	= [];
	while (<$fh>) {
		# Strip CRLF and LF
		s/[\r\n]+$//;
		# Line format: <name>:<val1>:<val2>:...
		my ($name, @values) = split(/:/, $_);
#		my $vector = vector(@values);
		my $vector = pdl( @values );
		# Store name and values
		push(@$names, 	$name);
		push(@$vectors,	$vector);
	}
	close($fh);
	# Return arrayrefs of names and vectors
	return ($names, $vectors);
}
############################################################
### Main Code
# Read command-line arguments
my ($DBfile, $matchfile, $outfile) = @ARGV;
die "Usage: $0 <DBfile> <matchfile> <outfile>\n" unless -f $DBfile && -f $matchfile && defined $outfile;
# Read inputs
my ($DBNames,		$DBVectors)		= readFile($DBfile);
my ($matchNames,	$matchVectors)	= readFile($matchfile);
open(my $fhOut, '>', $outfile);
# Write headers to results file
print $fhOut join(',', @$DBNames), "\n";
print $fhOut join(',', @$matchNames), "\n";
my $dbLen = @$DBVectors;
foreach my $mv (@$matchVectors) {
	my %zrec = ();
	foreach my $i (0 .. ($dbLen-1)) {
		my $dv = $DBVectors->[$i];
		# Calculate correlation of matchVector vs this DBvector
		my $zc = $mv->corr($dv);
		$zrec{$zc} = $i;
	}
	my @zsort = (sort {$b <=> $a} keys %zrec);
	my @zfull = ();
	foreach my $z (@zsort) {
		my $str = sprintf("%.8f:%d",$z,$zrec{$z});
		push(@zfull,$str);
	}
	my $zstr = join(', ',@zfull);
	print $fhOut "$zstr\n";
}
close($fhOut);
