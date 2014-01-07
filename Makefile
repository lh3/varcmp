.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

all:varcmp.pdf

varcmp.pdf:varcmp.tex varcmp.bib qst2/qroc-CHM1-bt2.pdf qst2/qroc-CHM1-mem.pdf qst2/qroc-NA12878-bt2.pdf qst2/qroc-NA12878-bwa.pdf qst2/qroc-NA12878-mem.pdf qst2/qroc-realn-bqsr.pdf tree.pdf
		pdflatex varcmp; bibtex varcmp; pdflatex varcmp; pdflatex varcmp;

qst2/qroc-CHM1-bt2.eps qst2/qroc-CHM1-mem.eps qst2/qroc-NA12878-bt2.eps qst2/qroc-NA12878-bwa.eps qst2/qroc-NA12878-mem.eps qst2/qroc-realn-bqsr.eps:
		(cd qst2; gnuplot plot.gp)

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ varcmp.pdf alnroc-?e.{eps,pdf}
