current_analysis: analysis/results/fpkm_sum Website
	@echo "Nothing deeper yet"

analysis/results/fpkm_sum: analysis/summary.tsv | analysis/results
	@echo "All genes should have approximately the same sum of FPKMs"
	python -c "import pandas as pd; print pd.read_table('analysis/summary.tsv',index_col=0).sum(axis=0)" \
		| tee $@
	
analysis/results/complexity:
	@echo "Temporarily deprecated"


analysis/results:
	mkdir analysis/results

Website: analysis/summary.tsv
	cp analysis/summary.tsv Website/genes.cuff
	python Website/draw_to_gene.py Website/genes.cuff
	echo $(MELGFF) > Website/versions.txt
	echo $(MELFASTA) >> Website/versions.txt

