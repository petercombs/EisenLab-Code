current-analysis: analysis/results/fpkm_sum \
	analysis/results/complexity
	@echo "Nothing deeper yet"

analysis/results/fpkm_sum: analysis/summary.tsv | analysis/results
	@echo "All genes should have approximately the same sum of FPKMs"
	python -c "import pandas as pd; print pd.read_table('analysis/summary.tsv',index_col=0).sum(axis=0)" \
		| tee $@
	
analysis/results/complexity: analysis/summary.tsv | analysis/results 
	python CheckCoverage.py \
		Reference/mel_good.gtf \
		analysis/*/accepted_hits_sorted.bam \
		| tee $@_mel


analysis/results:
	mkdir analysis/results

Reference/peaks-25-dm2/uptodate : | Reference/peaks-25-dm2
	python StealBindingBeds.py
	date > Reference/peaks-25-dm2/uptodate

Reference/peaks-25-dm2:
	mkdir $@

Reference/peaks-25-dm3/uptodate : Reference/peaks-25-dm2/uptodate | Reference/peaks-25-dm3 Reference/unmapped
	parallel "liftOver {} \
				       prereqs/dm2ToDm3.over.chain.gz \
	                   {//}/../peaks-25-dm3/{/} {//}/../unmapped/{/}" \
						   ::: Reference/peaks-25-dm2/*.bed
	date > Reference/peaks-25-dm2/uptodate

Reference/peaks-25-dm3:
	mkdir $@

Reference/unmapped:
	mkdir $@


BDTNPVER=2.1
prereqs/current_bdtnp:
	wget -O prereqs/current_bdtnp.tgz http://bdtnp.lbl.gov/Fly-Net/archives/chipper/BDTNP_in_vivo_binding_Release.$(BDTNPVER).tar.gz
	tar -xzf prereqs/current_bdtnp.tgz
	cp -r prereqs/BDTNP_in_vivo_binding_Release.$(BDTNPVER)/Supplemental_Tables/ $@


Reference/bcd_peaks: prereqs/current_bdtnp
	for TF in bcd cad da dl gt hb hkb kni kr mad med run shn slp1 sna tll twi z; do \
		cat prereqs/current_bdtnp/$${TF}_[1-3]_* \
			| awk 'NR > 1 {print $$2":"$$6}' \
			> /tmp/$${TF}_peaks; \
		echo 'OldPeak	NewPeak' > Reference/$${TF}_peaks; \
		curl -L  'http://flybase.org/cgi-bin/coord_converter.html' \
			-H 'Referer: http://flybase.org/static_pages/downloads/COORD.html' \
			-F species=dmel \
			-F inr=4\
			-F outr=6\
			-F saveas=File\
			-F ids="" \
			-F idfile="@/tmp/$${TF}_peaks" \
			-F .submit=Go \
			| grep -v '?' \
			>> Reference/$${TF}_peaks; \
	done


Reference/zld_peaks: prereqs/journal.pgen.1002266.s005.xls
	perl -pe 's//\n/g' $< \
		| awk 'NR > 2 {print $$2":"$$5}' \
		> /tmp/zld_peaks
	echo 'OldPeak	NewPeak' > $@
	curl -L  'http://flybase.org/cgi-bin/coord_converter.html' \
		-H 'Referer: http://flybase.org/static_pages/downloads/COORD.html' \
		-F species=dmel \
		-F inr=5\
		-F outr=6\
		-F saveas=File\
		-F ids="" \
		-F idfile="@/tmp/zld_peaks" \
		-F .submit=Go \
		| grep -v '?' \
		>> $@

Reference/tss: $(MELGTF)
	cat $< \
		| python FindTSSs.py \
		| rev \
		| uniq -f 1 \
		| rev \
		> $@

