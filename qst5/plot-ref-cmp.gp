set t po eps enh co so 'Helvetica,20'

set xran [0:600]
set out 'ref-cmp-snp.eps'
set xlab "#CHM1 heterozygous SNPs ({/Symbol \264}10^3)"
set ylab "#NA12878 - #CHM1 hets ({/Symbol \264}10^6)"
set pointsize 1.5
plot 'ref-cmp.hist' u ($2*1e-3):($3*1e-6) w p ls 7 notitle, '' u ($2*1e-3):($3*1e-6):1 w labels offset .8,.5 left notitle

set out 'ref-cmp-indel.eps'
set xran [30:*]
set yran [*:425]
set xlab "#CHM1 heterozygous INDELs ({/Symbol \264}10^3)"
set ylab "#NA12878 - #CHM1 hets ({/Symbol \264}10^3)"
plot 'ref-cmp.hist' u ($4*1e-3):($5*1e-3) w p ls 7 notitle, '' u ($4*1e-3):($5*1e-3):1 w labels offset -.8,.5 right notitle
