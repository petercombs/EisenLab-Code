# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis-tophat

# Reference FASTA and GFF files from FlyBase and SGD
MELRELEASE = r5.55_FB2014_01
VIRRELEASE = r1.2_FB2012_01
MELVERSION = $(word 1, $(subst _FB, ,$(MELRELEASE)))
VIRVERSION = $(word 1, $(subst _FB, ,$(VIRRELEASE)))
MELDATE = $(word 2, $(subst _FB, ,$(MELRELEASE)))

MELFASTA = prereqs/dmel-all-chromosome-$(MELVERSION).fasta
VIRFASTA = prereqs/dvir-all-chromosome-$(VIRVERSION).fasta

REFDIR = Reference

MELFASTA2= $(REFDIR)/dmel_prepend.fasta
VIRFASTA2= $(REFDIR)/dvir_prepend.fasta

ORTHOLOGS = prereqs/gene_orthologs_fb_$(MELDATE).tsv

CERFASTA = prereqs/S288C_reference_sequence_R64-1-1_20110203.fsa
CERFASTA2= $(REFDIR)/scer_prepend.fasta

MELGFF   = prereqs/dmel-all-$(MELVERSION).gff
MELGTF   = $(REFDIR)/mel_good.gtf
VIRGFF   = prereqs/dvir-all-$(VIRVERSION).gff
VIRGTF   = $(REFDIR)/vir_good.gtf
CERGFF   = prereqs/saccharomyces_cerevisiae_R64-1-1_20110208.gff
MELVIRGTF= $(REFDIR)/melvir.gtf
MELVIRFASTA=$(REFDIR)/melvir.fa


all : $(ANALYSIS_DIR)/summary.tsv current-analysis

# Read the per-project make-file
include config.make
include analyze.make

.SECONDARY: 

$(ANALYSIS_DIR) :
	mkdir $(ANALYSIS_DIR)

$(ANALYSIS_DIR)/summary.tsv : MakeSummaryTable.py $(FPKMS) $(RUNCONFIG) | $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py $(ANALYSIS_DIR) 

$(ANALYSIS_DIR)/summary_in_all.tsv : MakeSummaryTable.py $(FPKMS) $(RUNCONFIG) | $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py \
		--in-subdirectory all \
		$(ANALYSIS_DIR) 

$(ANALYSIS_DIR)/summary_in_subset.tsv : $(ANALYSIS_DIR)/subset_count MakeSummaryTable.py $(FPKMS_SUBSET) $(RUNCONFIG) | $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py \
		--filename subset_genes.fpkm_tracking \
		--in-subdirectory subset \
		$(ANALYSIS_DIR) 

%/genes.fpkm_tracking : %/assigned_dmelR.bam $(MELGTF) $(MELFASTA2)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(@D) -u \
		--frag-bias-correct $(MELFASTA2) -G $(MELGTF) $<

%/all/genes.fpkm_tracking : %/accepted_hits_sorted.bam $(MELVIRGTF) $(MELVIRFASTA)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(@D) -u \
		--frag-bias-correct $(MELVIRFASTA) -G $(MELVIRGTF) $<

%/subset/subset_genes.fpkm_tracking : $(ANALYSIS_DIR)/subset_count %/subset/accepted_hits_sorted.bam $(MELVIRGTF) $(MELVIRFASTA)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(@D) -u \
		--frag-bias-correct $(MELVIRFASTA) -G $(MELVIRGTF) \
		$(@D)/accepted_hits_sorted.bam
	mv $(@D)/genes.fpkm_tracking $@

$(ANALYSIS_DIR)/subset_count: $(FPKMS)
	python SubSample.py | tee $@

%/accepted_hits_sorted.bam: %/accepted_hits.bam
	samtools sort $< $(@D)/accepted_hits_sorted
	samtools index $@

%/assigned_dmelR.bam : %/accepted_hits.bam AssignReads2.py
	samtools view -H $< \
		| grep -Pv 'SN:(?!dmel)' \
		> $(@D)/mel_only.header.sam
	python AssignReads2.py $(@D)/accepted_hits.bam
	samtools sort $(@D)/assigned_dmel.bam \
		$(@D)/assigned_dmel_sorted
	samtools reheader $(@D)/mel_only.header.sam \
		$(@D)/assigned_dmel_sorted.bam > $@
	rm $(@D)/assigned_dmel_sorted.bam
	samtools index $@


$(MELGTF): $(MELGFF) | $(REFDIR)
	gffread $< -E -T -o- | \
		awk '{print "dmel_"$$0}' | \
		grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' > \
		$@

$(VIRGTF): $(VIRGFF) $(ORTHOLOGS) | $(REFDIR)
	gffread $< -E -T -o- \
		| awk '{print "dvir_"$$0}' \
		| grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' \
		| grep 'FBgn' \
		| python FilterOrthologs.py $(ORTHOLOGS) \
		> $@

$(MELFASTA): | $(REFDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/fasta/dmel-all-chromosome-$(MELVERSION).fasta.gz
	gunzip $@.gz

$(MELGFF): | $(REFDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/gff/dmel-all-$(MELVERSION).gff.gz
	gunzip $@.gz

$(VIRFASTA): | $(REFDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_$(VIRRELEASE)/fasta/dvir-all-chromosome-$(VIRVERSION).fasta.gz
	gunzip $@.gz

$(VIRGFF): | $(REFDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_virilis/dvir_$(VIRRELEASE)/gff/dvir-all-$(VIRVERSION).gff.gz
	gunzip $@.gz

$(MELFASTA2): $(MELFASTA)| $(REFDIR)
	perl -pe 's/>/>dmel_/' $(MELFASTA) > $@

$(VIRFASTA2): $(VIRFASTA)| $(REFDIR)
	perl -pe 's/>/>dvir_/' $(VIRFASTA) > $@

$(CERFASTA2): $(CERFASTA)| $(REFDIR)
	perl -pe 's/>/>scer_/' $(CERFASTA) > $@

$(MELVIRFASTA): $(MELFASTA2) $(VIRFASTA2)| $(REFDIR)
	cat $(MELFASTA2) $(VIRFASTA2) > $@


Reference/DmelScer/Genome : | $(MELFASTA2) $(CERFASTA2)  $(MELGTF) Reference/DmelScer $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir Reference/DmelScer \
		--genomeFastaFiles $(MELFASTA2) $(CERFASTA2) \
		--sjdbGTFfile $(MELGTF)

$(MELVIRGTF): $(MELGTF) $(VIRGTF) | $(REFDIR)
	cat $^ > $@

Reference/DmelScer: | $(REFDIR)
	mkdir $@

Reference/DmelDvir/transcriptome : |  Reference/DmelDvir
	tophat --GTF $(MELVIRGTF) \
		--transcriptome-index $@ \
		$(REFDIR)/DmelDvir
	touch $@


$(REFDIR)/DmelDvir: | $(REFDIR)
	mkdir $@

$(REFDIR)/DmelDvir/Genome : $(MELVIRGTF) |  Reference/DmelDvir $(MELFASTA2) $(VIRFASTA2)  $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir Reference/DmelDvir \
		--genomeFastaFiles $(MELFASTA2) $(VIRFASTA2) \
		--sjdbGTFfile $(MELVIRGTF)

$(ORTHOLOGS) :
	wget -O $@.gz -i ftp.flybase.org/releases/FB$(MELDATE)/precomputed_files/genes/gene_orthologs_fb_$(MELDATE).tsv.gz
	gunzip $@.gz

$(REFDIR) :
	mkdir $@
Reference/DmelDvir:
	bowtie2-build $(MELVIRFASTA) $@
	mkdir $@
