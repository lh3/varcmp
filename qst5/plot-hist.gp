set t po eps co so enh "Helvetica,20"

set style line 1 lc rgb "#555555" lw 1;
set style line 2 lc rgb "#ffffff" lw 1;
set style line 3 lc rgb "#777777" lw 1;
set style line 4 lc rgb "#bbbbbb" lw 1;
set style line 5 lc rgb "#999999" lw 1;

set size 1, 2



set out "hist-snp.eps"

set key top left
set multiplot layout 2,1 offset 0, .04

set style histogram rowstacked
set boxwidth 0.8 relative
set style data histograms
#set style fill solid 1.0 border lt -1
set style fill pattern 7 border lt -1
set ylab "#CHM1 heterozygous SNPs ({/Symbol \264}10^3)"

unset xtics
set size 1, 1
set origin 0, 1.04
set tmargin 3
set bmargin 0
plot \
	"hist.dat" u ($2*1e-3):xtic(1) t 'Filtered by LC' ls 1 fill pat 7, \
	"" u ($4*1e-3) t 'Filtered by DP' ls 2 fill pat 8, \
	"" u ($6*1e-3) t 'Filtered by misc' ls 3 fill pat 9, \
	"" u ($8*1e-3) t 'Remained' ls 4 fill pat 10


set ylab "#NA12878 heterozygous SNPs ({/Symbol \264}10^6)"
set xtics rotate by 315 nomirror
set size 1, 1
set tmargin 0
set bmargin 3
set yran [1.8:2.75]
plot \
	"hist.dat" u ($16*1e-6):xtic(1) t 'Remained' ls 4 fill pat 10, \
	"" u ($14*1e-6) t 'Filtered by misc' ls 3 fill pat 9, \
	"" u ($12*1e-6) t 'Filtered by DP' ls 2 fill pat 8, \
	"" u ($10*1e-6) t 'Filtered by LC' ls 1 fill pat 7

unset multiplot



set out "hist-indel.eps"

set key top right
set yran [*:*]
set multiplot layout 2,1 offset 0, .04

set style histogram rowstacked
set boxwidth 0.8 relative
set style data histograms
#set style fill solid 1.0 border lt -1
set style fill pattern 7 border lt -1
set ylab "#CHM1 heterozygous InDels ({/Symbol \264}10^3)"

unset xtics
set size 1, 1
set origin 0, 1.04
set tmargin 3
set bmargin 0
plot \
	"hist.dat" u ($18*1e-3):xtic(1) t '{/Symbol \261}1bp, Filtered by LC' ls 5 fill pat 5, \
	"" u (($3-$18)*1e-3) t '>1bp, Filtered by LC' ls 1 fill pat 6, \
	"" u ($5*1e-3) t 'Filtered by DP' ls 2 fill pat 8, \
	"" u ($7*1e-3) t 'Filtered by misc' ls 3 fill pat 9, \
	"" u ($9*1e-3) t 'Remained' ls 4 fill pat 10


set ylab "#NA12878 heterozygous InDels ({/Symbol \264}10^3)"
set xtics rotate by 315 nomirror
set size 1, 1
set tmargin 0
set bmargin 3
set yran [0:725]
plot \
	"hist.dat" u ($17*1e-3):xtic(1) t 'Remained' ls 4 fill pat 10, \
	"" u ($15*1e-3) t 'Filtered by misc' ls 3 fill pat 9, \
	"" u ($13*1e-3) t 'Filtered by DP' ls 2 fill pat 8, \
	"" u (($11-$19)*1e-3) t '>1bp, Filtered by LC' ls 1 fill pat 6, \
	"" u ($19*1e-3) t '{/Symbol \261}1bp, Filtered by LC' ls 5 fill pat 5

unset multiplot
