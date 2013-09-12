GENETABLE = prereqs/gene_map_table_fb_2013_04.tsv

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

Website: analysis/summary.tsv | Website/imgs
	cp analysis/summary.tsv Website/genes.cuff
	echo `basename $(MELGFF)` > Website/versions.txt
	echo `basename $(MELFASTA)` >> Website/versions.txt
	echo `basename $(GENETABLE)` >> Website/versions.txt
	echo 'Made on' >> Website/versions.txt
	date >> Website/versions.txt
	cut -f -2 $(GENETABLE) > Website/gene_table.tsv
	python Website/draw_to_gene.py Website/genes.cuff

Website/imgs:
	mkdir Website/imgs
