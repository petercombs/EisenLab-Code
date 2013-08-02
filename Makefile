# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis

# Reference FASTA and GFF files from FlyBase and SGD
MELFASTA = Reference/dmel-all-chromosome-r5.52.fasta
CERFASTA = Reference/S288C_reference_sequence_R64-1-1_20110203.fsa
MELGFF   = Reference/dmel-all-r5.52.gff
MELGTF   = Reference/mel_good.gtf
CERGFF   = Reference/saccharomyces_cerevisiae_R64-1-1_20110208.gff

all : $(ANALYSIS_DIR)/summary.tsv
	make -f analysis.make

# Read the per-project make-file
include config.make



$(ANALYSIS_DIR) :
	mkdir $(ANALYSIS_DIR)

$(ANALYSIS_DIR)/summary.tsv : MakeSummaryTable.py $(FPKMS) 
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py $(ANALYSIS_DIR) 

$(ANALYSIS_DIR)/%/genes.fpkm_tracking : $(ANALYSIS_DIR)/%/assigned_dmel.bam $(MELGTF)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(ANALYSIS_DIR)/$* -u \
		--frag-bias-correct $(MELFASTA) -G $(MELGTF) $<


# $(ANALYSIS_DIR)/%/accepted_hits.bam : $(ANALYSIS_DIR)/%/Aligned.out.sam 
#	samtools view -bS  -o $@  $<
#	rm $(ANALYSIS_DIR)/$*/Aligned.out.sam
#	# This sam file is big, let's get rid of it

$(ANALYSIS_DIR)/%/assigned_dmel.bam : $(ANALYSIS_DIR)/%/accepted_hits.bam AssignReads2.py
	samtools view -H $< | grep -Pv 'SN:(?!dmel)' > $(ANALYSIS_DIR)/$*/mel_only.header.sam
	python AssignReads2.py $(ANALYSIS_DIR)/$*/accepted_hits.bam
	samtools sort $(ANALYSIS_DIR)/$*/assigned_dmel.bam \
		$(ANALYSIS_DIR)/$*/assigned_dmel_sorted
	samtools reheader $(ANALYSIS_DIR)/$*/mel_only.header.sam \
		$(ANALYSIS_DIR)/$*/assigned_dmel_sorted.bam > $@
	rm $(ANALYSIS_DIR)/$*/assigned_dmel_sorted.bam


$(MELGTF): $(MELGFF)
	gffread $< -E -T -o- | \
		awk '{print dmel_$$0}' | \
		grep -vP '(snoRNA|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' > \
		$@

Reference/DmelScer/Genome : Reference/scer.fa Reference/dmel.fa
	
