.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

all:varcmp.pdf

varcmp.pdf:varcmp.tex varcmp.bib qst/qroc-aln1.pdf qst/qroc-aln2.pdf qst/qroc-aln3.pdf qst/qroc-bwa1.pdf qst/qroc-bwa2.pdf tree.pdf
		pdflatex varcmp; bibtex varcmp; pdflatex varcmp; pdflatex varcmp;

qst/qroc-aln1.eps qst/qroc-aln2.eps qst/qroc-aln3.eps qst/qroc-bwa1.eps qst/qroc-bwa2.eps:qst/qroc.gp
		(cd qst; gnuplot qroc.gp)

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ varcmp.pdf alnroc-?e.{eps,pdf}
