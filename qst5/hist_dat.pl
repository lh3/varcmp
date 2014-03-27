#!/usr/bin/perl -w

use strict;
use warnings;

die("Usage: hist_LC_DP.pl <prefix>") if @ARGV < 1;

my $prefix = $ARGV[0];
my (@s, @t, $fh);

open($fh, "$prefix.q1C") || die;
$_ = <$fh>;
@t = split;
@s = split while (<$fh>);
close($fh);
my $n_LC1_indel = $s[21] + $s[23] - $t[21] - $t[23];
my $m_LC1_indel = $s[10] + $s[12] - $t[10] - $t[12];

open($fh, "$prefix.qLC") || die;
$_ = <$fh>;
@t = split;
@s = split while (<$fh>);
close($fh);
my $n_LC_rest_snp = $t[15] + $t[17];
my $n_LC_rest_indel = $t[21] + $t[23];
my $n_LC_flt_snp = $s[15] + $s[17] - $n_LC_rest_snp;
my $n_LC_flt_indel = $s[21] + $s[23] - $n_LC_rest_indel;
my $m_LC_rest_snp = $t[4] + $t[6];
my $m_LC_rest_indel = $t[10] + $t[12];
my $m_LC_flt_snp = $s[4] + $s[6] - $m_LC_rest_snp;
my $m_LC_flt_indel = $s[10] + $s[12] - $m_LC_rest_indel;

open($fh, "$prefix.dQU") || die;
@s = split while (<$fh>);
my $n_DP_rest_snp = $s[15] + $s[17];
my $n_DP_rest_indel = $s[21] + $s[23];
my $n_DP_flt_snp = $n_LC_rest_snp - $n_DP_rest_snp;
my $n_DP_flt_indel = $n_LC_rest_indel - $n_DP_rest_indel;
my $m_DP_rest_snp = $s[4] + $s[6];
my $m_DP_rest_indel = $s[10] + $s[12];
my $m_DP_flt_snp = $m_LC_rest_snp - $m_DP_rest_snp;
my $m_DP_flt_indel = $m_LC_rest_indel - $m_DP_rest_indel;
close($fh);

open($fh, "$prefix.mQU") || die;
@s = split while (<$fh>);
my $n_misc_flt_snp = $n_DP_rest_snp - $s[15] - $s[17];
my $n_misc_flt_indel = $n_DP_rest_indel - $s[21] - $s[23];
my $m_misc_flt_snp = $m_DP_rest_snp - $s[4] - $s[6];
my $m_misc_flt_indel = $m_DP_rest_indel - $s[10] - $s[12];
close($fh);

print join("\t", $prefix,
			$n_LC_flt_snp, $n_LC_flt_indel, $n_DP_flt_snp, $n_DP_flt_indel, $n_misc_flt_snp, $n_misc_flt_indel, $n_DP_rest_snp - $n_misc_flt_snp, $n_DP_rest_indel - $n_misc_flt_indel,
			$m_LC_flt_snp, $m_LC_flt_indel, $m_DP_flt_snp, $m_DP_flt_indel, $m_misc_flt_snp, $m_misc_flt_indel, $m_DP_rest_snp - $m_misc_flt_snp, $m_DP_rest_indel - $m_misc_flt_indel,
			$n_LC1_indel, $m_LC1_indel), "\n";
