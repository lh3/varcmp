#!/usr/bin/perl

use strict;
use warnings;

my (%h, %cnt);
while (<>) {
	next if /realn/;
	my @t = split;
	$t[0] =~ s/\./:/;
	$t[1] =~ s/\./:/;
	$h{$t[0]}{$t[1]} = [$t[4]/$t[2], $t[6]];
	$cnt{$t[0]} = $t[2];
}

my @a = sort keys(%h);
for my $a (@a) {
	print " & $a";
}
print " \\\\\n";
for my $a (@a) {
	print $a;
	my ($min, $max) = (1, 0);
	for my $b (@a) {
		next if $a eq $b;
		$min = $min < $h{$a}{$b}[0]? $min : $h{$a}{$b}[0];
		$max = $max > $h{$a}{$b}[0]? $max : $h{$a}{$b}[0];
	}
	for my $b (@a) {
		if ($a eq $b) {
			print " & {\\bf $cnt{$a}}";
		} else {
			if (int($h{$a}{$b}[0]*1000 + .499) >= int($max*1000 + .499) - 5) {
				printf(" & {\\bf %.3f}", $h{$a}{$b}[0]);
			} else {
				printf(" & %.3f", $h{$a}{$b}[0]);
			}
			printf("/%.1f", $h{$a}{$b}[1]);
		}
	}
	print " \\\\\n";
}
