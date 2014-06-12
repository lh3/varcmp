#!/usr/bin/perl -w

use strict;
use warnings;

while (<>) {
	next if /qAB=(-?\d+)/ && $1 < 30;
	next if /qFS=(-?\d+)/ && $1 < -20;
	next if /qDS=(-?\d+)/ && $1 < 1;
	print;
}
