# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis

# Reference FASTA and GFF files from FlyBase and SGD
MELRELEASE = r6.01_FB2014_04
MELMAJORVERSION = $(word 1, $(subst ., , $(MELRELEASE)))
MELVERSION = $(word 1, $(subst _FB, ,$(MELRELEASE)))
MELDATE = $(word 2, $(subst _FB, ,$(MELRELEASE)))

MELFASTA = prereqs/dmel-all-chromosome-$(MELVERSION).fasta

REFDIR = Reference

MELFASTA2= $(REFDIR)/dmel_prepend.fasta

ORTHOLOGS = prereqs/gene_orthologs_fb_$(MELDATE).tsv

MELGFF   = prereqs/dmel-all-$(MELVERSION).gff
MELGTF   = $(REFDIR)/mel_good.gtf

GENEMAPTABLE = gene_map_table_fb_$(MELDATE).tsv


all : $(ANALYSIS_DIR)/summary.tsv

genomes: Reference/Dmel/Genome
	echo "Genomes Made"


# Read the per-project make-file
include config.make
include analyze.make

.SECONDARY:

$(ANALYSIS_DIR) :
	mkdir $(ANALYSIS_DIR)

$(ANALYSIS_DIR)/summary.tsv : MakeSummaryTable.py $(FPKMS) $(RUNCONFIG) Makefile | $(ANALYSIS_DIR)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py \
       --params $(RUNCONFIG) \
	   --strip-low-reads 1000000 \
	   --mapped-bamfile accepted_hits_sorted.bam \
		$(ANALYSIS_DIR)

%/genes.fpkm_tracking : %/accepted_hits_sorted.bam $(MELGTF) $(MELFASTA2)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(@D) -u \
		--frag-bias-correct $(MELFASTA2) -G $(MELGTF) $<

%/accepted_hits_sorted.bam: %/accepted_hits.bam
	samtools sort $< $(@D)/accepted_hits_sorted
	samtools index $@

$(MELGTF): $(MELGFF) | $(REFDIR)
	gffread $< -E -T -o- | \
		awk '{print "dmel_"$$0}' | \
		grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' > \
		$@

$(MELFASTA): $(REFDIR)/$(MELMAJORVERSION) | $(REFDIR)
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/fasta/dmel-all-chromosome-$(MELVERSION).fasta.gz
	gunzip $@.gz

$(MELGFF): $(REFDIR)/$(MELVERSION) | $(REFDIR) 
	wget -O $@.gz ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_$(MELRELEASE)/gff/dmel-all-$(MELVERSION).gff.gz
	gunzip $@.gz

$(REFDIR)/Dmel/transcriptome : | $(REFDIR)/Dmel
	tophat --GTF $(MELGTF) \
		--transcriptome-index $@ \
		$(REFDIR)/Dmel
	touch $@


$(REFDIR)/Dmel/Genome : $(MELGTF) |  $(REFDIR)/Dmel $(MELFASTA2) $(REFDIR)
	STAR --runMode genomeGenerate --genomeDir $(REFDIR)/Dmel \
		--genomeFastaFiles $(MELFASTA2) \
		--sjdbGTFfile $(MELGTF)

$(ORTHOLOGS) :
	wget -O $@.gz -i ftp.flybase.org/releases/FB$(MELDATE)/precomputed_files/genes/gene_orthologs_fb_$(MELDATE).tsv.gz
	gunzip $@.gz

$(REFDIR) :
	mkdir $@

$(REFDIR)/Dmel:
	bowtie2-build --offrate 1 $(MELFASTA2) $@
	mkdir $@

$(MELFASTA2): $(MELFASTA) | $(REFDIR)
	perl -pe 's/>/>dmel_/' $(MELFASTA) > $@

$(GENEMAPTABLE):
	wget ftp://ftp.flybase.net/releases/$(MELDATE)/precomputed_files/genes/$(GENEMAPTABLE).gz \
		-O $(GENEMAPTABLE).gz
	gunzip $(GENEMAPTABLE).gz

$(REFDIR)/$(MELVERSION):
	touch $@

$(REFDIR)/$(MELMAJORVERSION):
	touch $@
