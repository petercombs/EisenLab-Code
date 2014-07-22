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

Reference/unmapped:
	mkdir $@
