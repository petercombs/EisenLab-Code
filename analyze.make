current-analysis: analysis/results/fpkm_sum analysis/results/fold_change.log analysis/results/mel_summary.log
	@echo "Nothing deeper yet"

analysis/results/fpkm_sum: analysis/summary.tsv | analysis/results
	@echo "All genes should have approximately the same sum of FPKMs"
	python -c "import pandas as pd; print pd.read_table('analysis/summary.tsv',index_col=0).sum(axis=0)" \
		| tee $@
	
analysis/results/complexity:
	@echo "Temporarily deprecated"

analysis/results/fold_change.log: analysis/summary_in_all.tsv
	python AnalyzeFoldChanges.py $^ > $@

analysis/results/mel_summary.log: analysis/summary.tsv
	python AnalyzeMelConsistency.py $^ > $@

analysis/results:
	mkdir analysis/results

Reference/unmapped:
	mkdir $@
