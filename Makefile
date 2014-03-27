.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

all:varcmp.pdf

varcmp.pdf:varcmp.tex varcmp.bib qst5/hist-snp.pdf qst5/hist-indel.pdf qst5/roc1-co4.pdf
		pdflatex varcmp; bibtex varcmp; pdflatex varcmp; pdflatex varcmp;

qst5/hist-snp.eps qst5/hist-indel.eps:qst5/plot-hist.gp
		(cd qst5; gnuplot plot-hist.gp)

qst5/roc1-co4.eps:qst5/plot-co-w-DP4.gp
		(cd qst5; gnuplot plot-co-w-DP4.gp)

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ varcmp.pdf alnroc-?e.{eps,pdf} qst5/hist-snp.* qst5/hist-indel.* qst5/roc1-co4.*
