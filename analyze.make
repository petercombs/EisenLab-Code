current-analysis: analysis/results/fpkm_sum \
	analysis/results/fold_change.log \
	analysis/results/mel_summary.log \
	analysis/results/complexity
	@echo "Nothing deeper yet"

.PHONY: analysis/results/fpkm_sum analysis/results/complexity \
	analysis/results/fold_change.log analysis/results/mel_summary.log

analysis/results/fpkm_sum: analysis/summary.tsv | analysis/results
	@echo "All genes should have approximately the same sum of FPKMs"
	python -c "import pandas as pd; print pd.read_table('analysis/summary.tsv',index_col=0).sum(axis=0)" \
		| tee $@
	
analysis/results/complexity: | analysis/results
	python CheckCoverage.py \
		Reference/melvir.gtf \
		analysis/*/accepted_hits_sorted.bam \
		| tee $@
	python CheckCoverage.py \
		Reference/mel_good.gtf \
		analysis/*/accepted_hits_sorted.bam \
		| tee $@_mel
	python CheckCoverage.py \
		Reference/vir_good.gtf \
		analysis/*/accepted_hits_sorted.bam \
		| tee $@_vir


analysis/results/fold_change.log: analysis/summary_in_all.tsv | analysis/results
	python AnalyzeFoldChanges.py $^ > $@

analysis/results/mel_summary.log: analysis/summary.tsv | analysis/results
	python AnalyzeMelConsistency.py $^ > $@

analysis/results:
	mkdir analysis/results

Reference/unmapped:
	mkdir $@
