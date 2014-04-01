.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:varcmp.pdf varcmp.ps

varcmp.pdf:varcmp.tex varcmp.bib hist-snp.pdf hist-indel.pdf roc1-co4.pdf
		pdflatex varcmp; bibtex varcmp; pdflatex varcmp; pdflatex varcmp;

varcmp.dvi:varcmp.tex varcmp.bib hist-snp.eps hist-indel.eps roc1-co4.eps indel-exam.eps venn2/ref.eps \
	venn2/CHM1-indel.eps venn2/CHM1-snp.eps venn2/NA12878-LC-indel.eps venn2/NA12878-LC-snp.eps venn2/NA12878-nLC-indel.eps venn2/NA12878-nLC-snp.eps
		latex varcmp; bibtex varcmp; latex varcmp; latex varcmp

varcmp.ps:varcmp.dvi
		dvips varcmp.dvi

qst5/hist-snp.eps qst5/hist-indel.eps:qst5/plot-hist.gp
		(cd qst5; gnuplot plot-hist.gp)

qst5/roc1-co4.eps:qst5/plot-co-w-DP4.gp
		(cd qst5; gnuplot plot-co-w-DP4.gp)

hist-snp.eps:qst5/hist-snp.eps
		ln -s $<

hist-indel.eps:qst5/hist-indel.eps
		ln -s $<

roc1-co4.eps:qst5/roc1-co4.eps
		ln -s $<

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ varcmp.pdf varcmp.ps varcmp.dvi indel-exam.eps qst5/hist-snp.* qst5/hist-indel.* qst5/roc1-co4.* hist-snp.* hist-indel.* roc1-co4.* venn2/*.eps
