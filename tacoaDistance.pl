#!/usr/bin/perl
use strict;
use warnings;
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
print "Running\n";
my ($DBfile, $matchfile, $outfile) = @ARGV;
print "DB: $DBfile\n";
die "Usage: $0 <DBfile> <matchfile> <outfile>\n" unless -f $DBfile && -f $matchfile && defined $outfile;
# Read inputs
my ($DBNames,		$DBVectors)		= readFile($DBfile);
my ($matchNames,	$matchVectors)	= readFile($matchfile);
print "Read files\n";
open(my $fhOut, '>', $outfile);
print "Opened output file $outfile\n";
# Write headers to results file
print $fhOut join(',', @$DBNames), "\n";
print $fhOut join(',', @$matchNames), "\n";
my $dbLen = @$DBVectors;
foreach my $mv (@$matchVectors) {
	my %drec = ();
	foreach my $i (0 .. ($dbLen-1)) {
		my $dv = $DBVectors->[$i];
		# Norm vectors
        my $nmv = norm($mv);
        my $ndv = norm($dv);
        # Compute distance
        my $c = inner($nmv, $ndv);
        inner($nmv, $ndv, $c);
        my $dc = 1.0 - $c;
        my $K = exp(  -1 * ($dc ** 2.0)) / 2.0;
		$drec{$K} = $i;
	}
	my @dsort = (sort {$b <=> $a} keys %drec);
	my @dfull = ();
	foreach my $d (@dsort) {
		my $str = sprintf("%.8f:%d",$d,$drec{$d});
		push(@dfull,$str);
	}
	my $dstr = join(', ',@dfull);
	print $fhOut "$dstr\n";
}
close($fhOut);
