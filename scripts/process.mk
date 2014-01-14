all: test $(prefix).flt.gz

test:
ifeq ($(prefix),)
	@echo "Usage: make -f process.mk prefix=blabla"
	@exit 1
endif

$(prefix).deovlp.gz:$(prefix).vcf.gz
	k8 bio8.js deovlp $< | grep ^[0-9] | bgzip > $@

$(prefix).gt.gz:$(prefix).deovlp.gz
	k8 bio8.js upd1gt $< | bgzip > $@

$(prefix).gt.cbs:$(prefix).gt.gz
	k8 filter.js -OQ3 $< | k8 cnv-cbs.js > $@  

$(prefix).gt.dp:$(prefix).gt.gz
	k8 filter.js -OQ3 $< | tabtk num -Qc3 > $@

$(prefix).gt.dup:$(prefix).gt.cbs
	k8 cbs2bed.js $< > $@

$(prefix).flt.gz:$(prefix).gt.gz $(prefix).gt.dup $(prefix).gt.dp
	k8 filter.js -d5 -D `perl -ane 'print $$F[8]+3*($$F[8]-$$F[6]),"\n"' $(prefix).gt.dp` -Q3 -c $(prefix).gt.dup -t mdust-v28-p1.bed.gz $< | bgzip > $@
