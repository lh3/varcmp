set t po eps "Helvetica,20"

set pointsize 1.4
set xtics 0.1
set key bot left
set xlab "Accumulative SNP counts (million)"
set ylab "Marginal ts/tv"
set xran [3:4]
set yran [0:2.4]

set out "qroc-aln1.eps"
plot "NA12878-pe.bt2.fb.fq0.qst" u ($2*1e-6):9 t "bt2:fb" ls 1 w lp, \
	"NA12878-pe.bwa.fb.fq0.qst" u ($2*1e-6):9 t "bwa:fb" ls 2 w lp, \
	"NA12878-pe.mem.fb.fq0.qst" u ($2*1e-6):9 t "mem:fb" ls 3 w lp, \
	"NA12878-pe.bt2.hc.fq0.qst" u ($2*1e-6):9 t "bt2:hc" ls 4 w lp, \
	"NA12878-pe.bwa.hc.fq0.qst" u ($2*1e-6):9 t "bwa:hc" ls 6 w lp, \
	"NA12878-pe.mem.hc.fq0.qst" u ($2*1e-6):9 t "mem:hc" ls 8 w lp

set out "qroc-aln2.eps"
plot "NA12878-pe.bt2.st1.fq0.qst" u ($2*1e-6):9 t "bt2:st" ls 1 w lp, \
	"NA12878-pe.bwa.st1.fq0.qst" u ($2*1e-6):9 t "bwa:st" ls 2 w lp, \
	"NA12878-pe.mem.st1.fq0.qst" u ($2*1e-6):9 t "mem:st" ls 3 w lp, \
	"NA12878-pe.bt2.ug.fq0.qst" u ($2*1e-6):9 t "bt2:ug" ls 4 w lp, \
	"NA12878-pe.bwa.ug.fq0.qst" u ($2*1e-6):9 t "bwa:ug" ls 6 w lp, \
	"NA12878-pe.mem.ug.fq0.qst" u ($2*1e-6):9 t "mem:ug" ls 8 w lp

set yran [1.85:2.15]
set ylab "Accumulative ts/tv"

set out "qroc-bwa1.eps"
plot "NA12878-pe.bwa.fb.fq0.qst" u ($2*1e-6):6 t "bwa:fb" ls 1 w lp, \
	"NA12878-pe.bwa-realn-pl.fb.fq0.qst" u ($2*1e-6):6 t "bwa-realn:fb" ls 2 w lp, \
	"NA12878-pe.bwa-realn-pl-bqsr.fb.fq0.qst" u ($2*1e-6):6 t "bwa-realn-bqsr:fb" ls 3 w lp, \
	"NA12878-pe.bwa.hc.fq0.qst" u ($2*1e-6):6 t "bwa:hc" ls 4 w lp, \
	"NA12878-pe.bwa-realn-pl.hc.fq0.qst" u ($2*1e-6):6 t "bwa-realn:hc" ls 6 w lp, \
	"NA12878-pe.bwa-realn-pl-bqsr.hc.fq0.qst" u ($2*1e-6):6 t "bwa-realn-bqsr:hc" ls 8 w lp

set out "qroc-bwa2.eps"
plot "NA12878-pe.bwa.st1.fq0.qst" u ($2*1e-6):6 t "bwa:st" ls 1 w lp, \
	"NA12878-pe.bwa-realn-pl.st1.fq0.qst" u ($2*1e-6):6 t "bwa-realn:st" ls 2 w lp, \
	"NA12878-pe.bwa-realn-pl-bqsr.st1.fq0.qst" u ($2*1e-6):6 t "bwa-realn-bqsr:st" ls 3 w lp, \
	"NA12878-pe.bwa.ug.fq0.qst" u ($2*1e-6):6 t "bwa:ug" ls 4 w lp, \
	"NA12878-pe.bwa-realn-pl.ug.fq0.qst" u ($2*1e-6):6 t "bwa-realn:ug" ls 6 w lp, \
	"NA12878-pe.bwa-realn-pl-bqsr.ug.fq0.qst" u ($2*1e-6):6 t "bwa-realn-bqsr:ug" ls 8 w lp
