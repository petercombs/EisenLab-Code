current_analysis: analysis/results/fpkm_sum
	@echo "Nothing deeper yet"

analysis/results/fpkm_sum: analysis/summary.tsv | analysis/results
	@echo "All genes should have approximately the same sum of FPKMs"
	python -c "import pandas as pd; print pd.read_table('analysis/summary.tsv',index_col=0).sum(axis=0)" \
		| tee $@
	
analysis/results/complexity:
	@echo "Temporarily deprecated"


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
