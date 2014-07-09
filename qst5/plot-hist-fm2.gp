set t po eps co so enh "Helvetica,20"

set style line 1 lc rgb "#555555" lw 1;
set style line 2 lc rgb "#ffffff" lw 1;
set style line 3 lc rgb "#777777" lw 1;
set style line 4 lc rgb "#bbbbbb" lw 1;
set style line 5 lc rgb "#999999" lw 1;

set style line 1 lt 1 lc rgb "#FF0000" lw 1;
set style line 2 lt 1 lc rgb "#00C000" lw 1;
set style line 3 lt 1 lc rgb "#0080FF" lw 1;
set style line 4 lt 1 lc rgb "#C000FF" lw 1;
set style line 5 lt 1 lc rgb "#00EEEE" lw 1;
set style line 6 lt 1 lc rgb "#FF80FF" lw 1;

set size 1, 2
set key invert


set out "hist-snp-fm2.eps"

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
	"hist-fm2.dat" u ($8*1e-3):xtic(1) t 'Remained' ls 4 fill pat 10, \
	"" u ($6*1e-3) t 'Filtered by misc' ls 2 fill pat 9, \
	"" u ($4*1e-3) t 'Filtered by DP' ls 3 fill pat 7, \
	"" u ($2*1e-3) t 'Filtered by LC' ls 1 fill pat 6


set ylab "#NA12878 - #CHM1 hets ({/Symbol \264}10^6)"
set xtics rotate by 315 nomirror
set size 1, 1
set tmargin 0
set bmargin 3
set key top right
set yran [1.8:2.37]
plot \
	"hist-fm2.dat" u (($16-$8)*1e-6):xtic(1) t 'Remained' ls 4 fill pat 10, \
	"" u (($14-$6)*1e-6) t 'Filtered by misc' ls 2 fill pat 9, \
	"" u (($12-$4)*1e-6) t 'Filtered by DP' ls 3 fill pat 7, \
	"" u (($10-$2)*1e-6) t 'Filtered by LC' ls 1 fill pat 6

unset multiplot



set out "hist-indel-fm2.eps"

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
	"hist-fm2.dat" u ($9*1e-3):xtic(1) t 'Remained' ls 4 fill pat 10, \
	"" u ($7*1e-3) t 'Filtered by misc' ls 2 fill pat 9, \
	"" u ($5*1e-3) t 'Filtered by DP' ls 3 fill pat 7, \
	"" u (($3-$18)*1e-3) t '>1bp, Filtered by LC' ls 1 fill pat 6, \
	"" u ($18*1e-3):xtic(1) t '{/Symbol \261}1bp, Filtered by LC' ls 6 fill pat 5


set ylab "#NA12878 - #CHM1 hets ({/Symbol \264}10^3)"
set xtics rotate by 315 nomirror
set size 1, 1
set tmargin 0
set bmargin 3
set yran [0:550]
plot \
	"hist-fm2.dat" u (($17-$9)*1e-3):xtic(1) t 'Remained' ls 4 fill pat 10, \
	"" u (($15-$7)*1e-3) t 'Filtered by misc' ls 2 fill pat 9, \
	"" u (($13-$3)*1e-3) t 'Filtered by DP' ls 3 fill pat 7, \
	"" u (($11-$19-($3-$18))*1e-3) t '>1bp, Filtered by LC' ls 1 fill pat 6, \
	"" u (($19-$18)*1e-3) t '{/Symbol \261}1bp, Filtered by LC' ls 6 fill pat 5

unset multiplot
