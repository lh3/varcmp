.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:varcmp.pdf varcmp.ps

varcmp.pdf:varcmp.tex varcmp.bib hist-snp.pdf hist-indel.pdf roc1-co4.pdf \
	ref.pdf CHM1-indel.pdf CHM1-snp.pdf NA12878-LC-indel.pdf NA12878-LC-snp.pdf NA12878-nLC-indel.pdf NA12878-nLC-snp.pdf
		pdflatex varcmp; bibtex varcmp; pdflatex varcmp; pdflatex varcmp;

varcmp.dvi:varcmp.tex varcmp.bib hist-snp.eps hist-indel.eps roc1-co4.eps indel-exam.eps \
	ref.eps CHM1-indel.eps CHM1-snp.eps NA12878-LC-indel.eps NA12878-LC-snp.eps NA12878-nLC-indel.eps NA12878-nLC-snp.eps
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

ref.pdf:
		ln -s venn2/ref.pdf

CHM1-indel.pdf:venn2/CHM1-indel.pdf
		ln -s $<

CHM1-snp.pdf:venn2/CHM1-snp.pdf
		ln -s $<

NA12878-LC-indel.pdf:venn2/NA12878-LC-indel.pdf
		ln -s $<

NA12878-LC-snp.pdf:venn2/NA12878-LC-snp.pdf
		ln -s $<

NA12878-nLC-indel.pdf:venn2/NA12878-nLC-indel.pdf
		ln -s $<

NA12878-nLC-snp.pdf:venn2/NA12878-nLC-snp.pdf
		ln -s $<

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ varcmp.pdf varcmp.ps varcmp.dvi indel-exam.eps qst5/hist-snp.* qst5/hist-indel.* qst5/roc1-co4.* hist-snp.* hist-indel.* roc1-co4.* \
			*.eps NA12878-*.pdf CHM1-*.pdf ref.pdf ref.eps
