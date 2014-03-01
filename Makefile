# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis

# Reference FASTA and GFF files from FlyBase and SGD
MELRELEASE = r5.55_FB2014_01
VIRRELEASE = r1.2_FB2012_01
MELVERSION = $(word 1, $(subst _FB, ,$(MELRELEASE)))
VIRVERSION = $(word 1, $(subst _FB, ,$(VIRRELEASE)))
MELDATE = $(word 2, $(subst _FB, ,$(MELRELEASE)))

MELFASTA = prereqs/dmel-all-chromosome-$(MELVERSION).fasta
VIRFASTA = prereqs/dvir-all-chromosome-$(VIRVERSION).fasta

MELFASTA2= Reference/dmel_prepend.fasta
VIRFASTA2= Reference/dvir_prepend.fasta

ORTHOLOGS = prereqs/gene_orthologs_fb_$(MELDATE).tsv

CERFASTA = prereqs/S288C_reference_sequence_R64-1-1_20110203.fsa
CERFASTA2= Reference/scer_prepend.fasta

MELGFF   = prereqs/dmel-all-$(MELVERSION).gff
MELGTF   = Reference/mel_good.gtf
VIRGFF   = prereqs/dvir-all-$(VIRVERSION).gff
VIRGTF   = Reference/vir_good.gtf
CERGFF   = prereqs/saccharomyces_cerevisiae_R64-1-1_20110208.gff
MELVIRGTF= Reference/melvir.gtf


all : $(ANALYSIS_DIR)/summary.tsv current-analysis

# Read the per-project make-file
include config.make
include analyze.make

.SECONDARY: 

$(ANALYSIS_DIR) :
	mkdir $(ANALYSIS_DIR)

$(ANALYSIS_DIR)/summary.tsv : MakeSummaryTable.py $(FPKMS) $(RUNCONFIG)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py $(ANALYSIS_DIR) 

$(ANALYSIS_DIR)/%/genes.fpkm_tracking : $(ANALYSIS_DIR)/%/assigned_dmelR.bam $(MELGTF) $(MELFASTA2)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(ANALYSIS_DIR)/$* -u \
		--frag-bias-correct $(MELFASTA2) -G $(MELGTF) $<


# $(ANALYSIS_DIR)/%/accepted_hits.bam : $(ANALYSIS_DIR)/%/Aligned.out.sam 
#	samtools view -bS  -o $@  $<
#	rm $(ANALYSIS_DIR)/$*/Aligned.out.sam
#	# This sam file is big, let's get rid of it

$(ANALYSIS_DIR)/%/assigned_dmelR.bam : $(ANALYSIS_DIR)/%/accepted_hits.bam AssignReads2.py
	samtools view -H $< \
		| grep -Pv 'SN:(?!dmel)' \
		> $(ANALYSIS_DIR)/$*/mel_only.header.sam
	python AssignReads2.py $(ANALYSIS_DIR)/$*/accepted_hits.bam
	samtools sort $(ANALYSIS_DIR)/$*/assigned_dmel.bam \
		$(ANALYSIS_DIR)/$*/assigned_dmel_sorted
	samtools reheader $(ANALYSIS_DIR)/$*/mel_only.header.sam \
		$(ANALYSIS_DIR)/$*/assigned_dmel_sorted.bam > $@
	rm $(ANALYSIS_DIR)/$*/assigned_dmel_sorted.bam
	samtools index $@


$(MELGTF): $(MELGFF)
	gffread $< -E -T -o- | \
		awk '{print "dmel_"$$0}' | \
		grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' > \
		$@

$(VIRGTF): $(VIRGFF) $(ORTHOLOGS)
	gffread $< -E -T -o- \
		| awk '{print "dvir_"$$0}' \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' \
		| grep 'FBgn' \
		| python FilterOrthologs.py $(ORTHOLOGS) \
		> $@

$(MELFASTA): 
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/fasta/dmel-all-chromosome-$(MELVERSION).fasta.gz
	gunzip $@.gz

$(MELGFF): 
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/gff/dmel-all-$(MELVERSION).gff.gz
	gunzip $@.gz

$(VIRFASTA): 
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_$(VIRRELEASE)/fasta/dvir-all-chromosome-$(VIRVERSION).fasta.gz
	gunzip $@.gz

$(VIRGFF): 
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_$(VIRRELEASE)/gff/dvir-all-$(VIRVERSION).gff.gz
	gunzip $@.gz

$(MELFASTA2): $(MELFASTA)
	perl -pe 's/>/>dmel_/' $(MELFASTA) > $@

$(VIRFASTA2): $(VIRFASTA)
	perl -pe 's/>/>dmel_/' $(VIRFASTA) > $@

$(CERFASTA2): $(CERFASTA)
	perl -pe 's/>/>scer_/' $(CERFASTA) > $@



Reference/DmelScer/Genome : | $(MELFASTA2) $(CERFASTA2)  $(MELGTF) Reference/DmelScer
	STAR --runMode genomeGenerate --genomeDir Reference/DmelScer \
		--genomeFastaFiles $(MELFASTA2) $(CERFASTA2) \
		--sjdbGTFfile $(MELGTF)

$(MELVIRGTF): $(MELGTF) $(VIRGTF)
	cat $^ > $@

Reference/DmelScer:
	mkdir $@


Reference/DmelDvir:
	mkdir $@

Reference/DmelDvir/Genome : $(MELVIRGTF) |  Reference/DmelDvir $(MELFASTA2) $(VIRFASTA2) 
	STAR --runMode genomeGenerate --genomeDir Reference/DmelDvir \
		--genomeFastaFiles $(MELFASTA2) $(VIRFASTA2) \
		--sjdbGTFfile $(MELVIRGTF)

$(ORTHOLOGS) :
	wget -O $@.gz -i ftp.flybase.org/releases/FB$(MELDATE)/precomputed_files/genes/gene_orthologs_fb_$(MELDATE).tsv.gz
	gunzip $@.gz

