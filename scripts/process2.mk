all:test $(dir)/CHM1.$(mid).flt.qual $(dir)/CHM1.$(mid).flt.qAB $(dir)/CHM1.$(mid).flt.qDP \
	$(dir)/CHM1.$(mid).flt.qDS $(dir)/CHM1.$(mid).flt.qFS $(dir)/CHM1.$(mid).flt.qDS $(dir)/CHM1.$(mid).flt.qHW \
	$(dir)/CHM1.$(mid).flt.qLC \
	$(dir)/CHM1.$(mid).flt.dAB $(dir)/CHM1.$(mid).flt.dQU $(dir)/CHM1.$(mid).flt.dDS \
	$(dir)/CHM1.$(mid).flt.dFS $(dir)/CHM1.$(mid).flt.dHW \
	$(dir)/CHM1.$(mid).flt.mQU \
	$(dir)/CHM1.$(mid).flt.q1C

dir='./'

test:
ifeq ($(mid),)
	@echo "Usage: make -f process2.mk dir=foo mid=bar"
	@exit 1
endif

$(dir)/CHM1.$(mid).flt.qual:
	k8 pair-qst1.js $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.qAB:
	k8 pair-qst1.js -t qAB $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.qDP:
	k8 pair-qst1.js -t qDP $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.qDS:
	k8 pair-qst1.js -t qDS $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.qFS:
	k8 pair-qst1.js -t qFS $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.qHW:
	k8 pair-qst1.js -t qHW $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@


$(dir)/CHM1.$(mid).flt.qLC:
	k8 pair-qst1.js -t qLC $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.q1C:
	k8 pair-qst1.js -1t qLC $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@


$(dir)/CHM1.$(mid).flt.dQU:
	k8 pair-qst1.js -D $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.dAB:
	k8 pair-qst1.js -Dt qAB $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.dDS:
	k8 pair-qst1.js -Dt qDS $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.dFS:
	k8 pair-qst1.js -Dt qFS $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.dHW:
	k8 pair-qst1.js -Dt qHW $(dir)/NA12878.$(mid).flt.gz $(dir)/CHM1.$(mid).flt.gz > $@

$(dir)/CHM1.$(mid).flt.mQU:
	bash -c 'k8 pair-qst1.js -D <(zcat $(dir)/NA12878.$(mid).flt.gz|perl misc-flt.pl) <(zcat $(dir)/CHM1.$(mid).flt.gz|perl misc-flt.pl) > $@'

