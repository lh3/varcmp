set t po eps  enh "Helvetica,18"

set style line 1 pt 1 lw 2;
set style line 2 pt 2 lw 2;
set style line 3 pt 4 lw 2;
set style line 4 pt 6 lw 2;
set style line 6 pt 8 lw 2;

set pointsize 1.4

set out "roc1-co4.eps"

set yran [1.85:2.15]
set xran [0:140]
set key top right

set size 0.86, 1.2

set multiplot layout 2,2 offset .12,.11

set xtics format ""
set ylab "#NA12878 hets - #CHM1 hets ({/Symbol \264}10^6)"
set ytics .1

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set title "fb"
set size .35,.5
set arrow from 33,2.09 to 41,2.05
plot 'CHM1.bt2.fb.flt.qual' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) t 'QU' w lp, \
	 'CHM1.bt2.fb.flt.qDP' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) t 'DP' ls 2 w lp, \
	 'CHM1.bt2.fb.flt.qAB' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) t 'AB' ls 3 w lp, \
	 'CHM1.bt2.fb.flt.qFS' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) t 'FS' ls 4 w lp, \
	 'CHM1.bt2.fb.flt.qDS' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) t 'DS' ls 6 w lp
set ytics format ""
unset arrow
unset ylab

set y2lab "bt2"
set title "ug"
set size .35,.5
set origin .47,.61
set arrow from 69,2.14 to 87,2.12
plot 'CHM1.bt2.ug.flt.qual' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 1 w lp, \
	 'CHM1.bt2.ug.flt.qDP' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 2 w lp, \
	 'CHM1.bt2.ug.flt.qAB' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 3 w lp, \
	 'CHM1.bt2.ug.flt.qFS' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 4 w lp
unset arrow
unset y2lab


set xlab "#heterozygous SNPs in CHM1 ({/Symbol \264}10^3)"
set xtics 20
set xtics format "%.0f"
unset title

set size .35,.5
set ytics .1
set arrow from 43,2.13 to 51,2.10
plot 'CHM1.mem.fb.flt.qual' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 1 w lp, \
	 'CHM1.mem.fb.flt.qDP' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 2 w lp, \
	 'CHM1.mem.fb.flt.qAB' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 3 w lp, \
	 'CHM1.mem.fb.flt.qFS' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 4 w lp, \
	 'CHM1.mem.fb.flt.qDS' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 6 w lp
unset arrow

set y2lab "mem"
set xtics format ""
unset xlab
set size .35,.5
set origin .47,.11
set arrow from 45,2.13 to 63,2.102
plot 'CHM1.mem.ug.flt.qual' u (($16+$18)*1e-3):(($5+$7-$16-$18)*1e-6) not ls 1 w lp, \
	 'CHM1.mem.ug.flt.qDP' u (($16+$18)*1e-3):(($5+$7-$16-$18)* 1e-6) not ls 2 w lp, \
	 'CHM1.mem.ug.flt.qAB' u (($16+$18)*1e-3):(($5+$7-$16-$18)* 1e-6) not ls 3 w lp, \
	 'CHM1.mem.ug.flt.qFS' u (($16+$18)*1e-3):(($5+$7-$16-$18)* 1e-6) not ls 4 w lp

unset multiplot
