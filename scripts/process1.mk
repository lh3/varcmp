all:test $(prefix).flt.gz $(prefix).flt.indel1.gz $(prefix).flt.indel2.gz $(prefix).flt.snp1.gz $(prefix).flt.snp2.gz

test:
ifeq ($(prefix),)
	@echo "Usage: make -f process1.mk prefix=blabla"
	@exit 1
endif

$(prefix).deovlp.gz:$(prefix).vcf.gz
	k8 bio8.js deovlp $< | egrep "^(chr)?[0-9]" | bgzip > $@

$(prefix).gt.gz:$(prefix).deovlp.gz
	k8 bio8.js upd1gt $< | bgzip > $@

$(prefix).flt.gz:$(prefix).gt.gz
	k8 vcf-extra-info.js -h 1000g.hwe-bad.bed -l LCR-hs37d5.bed.gz $< | bgzip > $@

$(prefix).flt.indel1.gz:$(prefix).flt.gz
	k8 vcf-extra-flt.js -L $< | grep "LC=0;" | k8 bio8.js qst1 -G | bgzip > $@

$(prefix).flt.indel2.gz:$(prefix).flt.gz
	k8 vcf-extra-flt.js -L $< | grep "LC=-[0-9]*;" | k8 bio8.js qst1 -G | bgzip > $@

$(prefix).flt.snp1.gz:$(prefix).flt.gz
	k8 vcf-extra-flt.js -L $< | grep "LC=0;" | k8 bio8.js qst1 -S | bgzip > $@

$(prefix).flt.snp2.gz:$(prefix).flt.gz
	k8 vcf-extra-flt.js -L $< | grep "LC=-[0-9]*;" | k8 bio8.js qst1 -S | bgzip > $@
