#!/usr/bin/env perl

## Quick-and-dirty way to prepend copyleft.

$#ARGV >= 1 || die "usage: prepend-notice.pl <notice> <file 1> [<file 2> ...]\n";

$copying_src = shift @ARGV;

print "## applying notice in $copying_src\n";

$copying = `cat $copying_src`;

foreach (@ARGV) {
    if (/^${copying_src}$/) {
	print STDERR "skipping $copying_src\n";
    }
    else {
	print "## processing $_\n";
	$orig = `cat $_`;
	open(FOUT, ">$_");
	print FOUT "% $_\n\n";
	print FOUT "$copying\n";
	print FOUT "$orig";
	close(FOUT);
    }
}
