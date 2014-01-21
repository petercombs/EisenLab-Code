GENETABLE = prereqs/gene_map_table_fb_2013_04.tsv

current-analysis: analysis/results/fpkm_sum Website 
	@echo "Nothing deeper yet"

analysis/results/fpkm_sum: analysis/summary.tsv | analysis/results
	@echo "All genes should have approximately the same sum of FPKMs"
	python -c "import pandas as pd; print pd.read_table('analysis/summary.tsv',index_col=0).sum(axis=0)" \
		| tee $@
	
analysis/results/complexity:
	@echo "Temporarily deprecated"


analysis/results:
	mkdir analysis/results

Website: analysis/summary.tsv Website/draw_to_gene.py | Website/imgs
	cp analysis/summary.tsv Website/genes.cuff
	echo `basename $(MELGFF)` > Website/versions.txt
	echo `basename $(MELFASTA)` >> Website/versions.txt
	echo `basename $(GENETABLE)` >> Website/versions.txt
	echo 'Made on' >> Website/versions.txt
	date >> Website/versions.txt
	cut -f -2 $(GENETABLE) > Website/gene_table.tsv
	python Website/draw_to_gene.py Website/genes.cuff

FlyBase/FlyBase.tgz: analysis/summary_flybase.tsv | FlyBase FlyBase/imgs
	cp analysis/summary_flybase.tsv FlyBase
	echo `basename $(MELGFF)` > FlyBase/versions.txt
	echo `basename $(MELFASTA)` >> FlyBase/versions.txt
	echo `basename $(GENETABLE)` >> FlyBase/versions.txt
	echo 'Made on' >> FlyBase/versions.txt
	date >> FlyBase/versions.txt
	rm -r FlyBase/imgs
	mkdir FlyBase/imgs
	python Website/draw_to_gene.py -W 550 -H 75 -w 175 FlyBase/summary_flybase.tsv
	tar -cvzf FlyBase/FlyBase.tgz \
		FlyBase/imgs FlyBase/versions.txt FlyBase/summary_flybase.tsv

analysis/summary_flybase.tsv:
	if [ -a analysis/summary.tsv]; \
		then mv analysis/summary.tsv analysis/tmptmptmp.tmp; fi;
	python MakeSummaryTable.py --params $(RUNCONFIG) --key tracking_id \
		$(ANALYSIS_DIR)
	mv analysis/summary.tsv analysis/summary_flybase.tsv
	if [ -a analysis/tmptmptmp.tmp]; \
		then mv analysis/tmptmptmp.tmp analysis/summary.tsv; fi;


FlyBase:
	mkdir FlyBase

FlyBase/imgs: | FlyBase
	mkdir FlyBase/imgs

Website/imgs:
	mkdir Website/imgs
