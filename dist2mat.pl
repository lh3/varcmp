#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts("j", \%opts);

my $is_jdist = defined($opts{j});
my %h;

while (<>) {
	s/NA12878-pe.//g;
	s/.flt.gz//g;
	s/-pl//g;
	my @t = split;
	$h{$t[0]}{$t[1]} = 1 - ($is_jdist? $t[4] / ($t[2] < $t[3]? $t[2] : $t[3]) : $t[4] / ($t[2] + $t[3] - $t[4]));
}

my @a = keys(%h);
print scalar(@a), "\n";
for my $a (@a) {
	print $a;
	for my $b (@a) {
		printf("\t%.6f", $h{$a}{$b});
	}
	print "\n";
}
