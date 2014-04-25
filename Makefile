# Configuration files for the experiment
RUNCONFIG  = Parameters/RunConfig.cfg
STARCONFIG = Parameters/STAR_params.in

# Other random variables
ANALYSIS_DIR = analysis-multi

# Reference FASTA and GFF files from FlyBase and SGD
MELFASTA = prereqs/dmel-all-chromosome-r5.53.fasta
MELFASTA2= Reference/dmel_prepend.fasta
CERFASTA = prereqs/S288C_reference_sequence_R64-1-1_20110203.fsa
CERFASTA2= Reference/scer_prepend.fasta
MELGFF   = prereqs/dmel-all-r5.53.gff
MELGTF   = Reference/mel_good.gtf
CERGFF   = prereqs/saccharomyces_cerevisiae_R64-1-1_20110208.gff


all : $(ANALYSIS_DIR)/summary.tsv 

genomes: Reference/DmelDwil/Genome Reference/DmelDvir/Genome Reference/DmelDper/Genome Reference/DmelDmoj/Genome
	echo "Genomes Made"


# Read the per-project make-file
include config.make
include analyze.make

ambigs: $(FPKMS_AMBIG)

.SECONDARY: 

$(ANALYSIS_DIR) :
	mkdir $(ANALYSIS_DIR)

$(ANALYSIS_DIR)/summary.tsv : MakeSummaryTable.py $(FPKMS) $(RUNCONFIG)
	@echo '============================='
	@echo 'Making summary table'
	@echo '============================='
	python MakeSummaryTable.py --params $(RUNCONFIG) $(ANALYSIS_DIR) 

$(ANALYSIS_DIR)/%/genes.fpkm_tracking : $(ANALYSIS_DIR)/%/assigned_dmelR.bam $(MELGTF) $(MELFASTA2)
	@echo '============================='
	@echo 'Calculating Abundances'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(@D) -u \
		--frag-bias-correct $(MELFASTA2) -G $(MELGTF) $<

$(ANALYSIS_DIR)/%/withambig/genes.fpkm_tracking : $(ANALYSIS_DIR)/%/dmel_ambig_merged.bam $(MELGTF) $(MELFASTA2)
	@echo '============================='
	@echo 'Calculating Abundances with Ambiguity'
	@echo '============================='
	cufflinks --num-threads 8 --output-dir $(@D) -u \
		--frag-bias-correct $(MELFASTA2) -G $(MELGTF) $<

# $(ANALYSIS_DIR)/%/accepted_hits.bam : $(ANALYSIS_DIR)/%/Aligned.out.sam 
#	samtools view -bS  -o $@  $<
#	rm $(ANALYSIS_DIR)/$*/Aligned.out.sam
#	# This sam file is big, let's get rid of it

$(ANALYSIS_DIR)/%/dmel_ambig_merged.bam : $(ANALYSIS_DIR)/%/assigned_dmelR.bam
	echo $(@D)
	samtools view $(@D)/ambiguous.bam \
		| grep 'dmel' \
		| cat $(@D)/mel_only.header.sam - \
		| samtools view -bS -o $(@D)/ambig_dmelR.bam -
	samtools sort $(@D)/ambig_dmelR.bam $(@D)/ambig_dmelRS
	samtools merge $@ $< $(@D)/ambig_dmelRS.bam
	rm $(@D)/ambig_dmelR.bam

$(ANALYSIS_DIR)/%/assigned_dmelR.bam : $(ANALYSIS_DIR)/%/accepted_hits.bam AssignReads2.py
	samtools view -H $< | grep -Pv 'SN:(?!dmel)' > $(@D)/mel_only.header.sam
	python AssignReads2.py $(@D)/accepted_hits.bam
	samtools sort $(@D)/assigned_dmel.bam \
		$(@D)/assigned_dmel_sorted
	samtools view $(@D)/assigned_dmel_sorted.bam \
		| cat $(@D)/mel_only.header.sam - \
		| samtools view -bS -o $@ -
	rm $(@D)/assigned_dmel_sorted.bam
	samtools index $@


$(MELGTF): $(MELGFF)
	gffread $< -E -T -o- | \
		awk '{print "dmel_"$$0}' | \
		grep -vP '(snoRNA|CR[0-9]{4}|Rp[ILS]|mir-|tRNA|unsRNA|snRNA|snmRNA|scaRNA|rRNA|RNA:|mt:)' > \
		$@

$(MELFASTA2): $(MELFASTA)
	perl -pe 's/>/>dmel_/' $(MELFASTA) > $@

$(CERFASTA2): $(CERFASTA)
	perl -pe 's/>/>scer_/' $(CERFASTA) > $@

Reference/DmelScer/Genome : | $(MELFASTA2) $(CERFASTA2)  $(MELGTF) Reference/DmelScer
	STAR --runMode genomeGenerate --genomeDir Reference/DmelScer \
		--genomeFastaFiles $(MELFASTA2) $(CERFASTA2) \
		--sjdbGTFfile $(MELGTF)

Reference/DmelDper/transcriptome : |  Reference/DmelDper
	tophat --GTF Reference/AAA/melper.gtf \
		--transcriptome-index $@ \
		Reference/AAA/melper
	touch $@

Reference/DmelDwil/transcriptome : |  Reference/DmelDwil
	tophat --GTF Reference/AAA/melwil.gtf \
		--transcriptome-index $@ \
		Reference/AAA/melwil
	touch $@

Reference/DmelDmoj/transcriptome : |  Reference/DmelDmoj
	tophat --GTF Reference/AAA/melmoj.gtf \
		--transcriptome-index $@ \
		Reference/AAA/melmoj
	touch $@

Reference/DmelDvir/transcriptome : |  Reference/DmelDvir
	tophat --GTF Reference/AAA/melvir.gtf \
		--transcriptome-index $@ \
		Reference/AAA/melvir
	touch $@

Reference/DmelDper/Genome : |  Reference/DmelDper
	STAR --runMode genomeGenerate --genomeDir Reference/DmelDper \
		--genomeFastaFiles Reference/AAA/melper.fa \
		--sjdbGTFfile Reference/AAA/melper.gtf

Reference/DmelDwil/Genome : |  Reference/DmelDwil
	STAR --runMode genomeGenerate --genomeDir Reference/DmelDwil \
		--genomeFastaFiles Reference/AAA/melwil.fa \
		--sjdbGTFfile Reference/AAA/melwil.gtf

Reference/DmelDvir/Genome : |  Reference/DmelDvir
	STAR --runMode genomeGenerate --genomeDir Reference/DmelDvir \
		--genomeFastaFiles Reference/AAA/melvir.fa \
		--sjdbGTFfile Reference/AAA/melvir.gtf

Reference/DmelDmoj/Genome : |  Reference/DmelDmoj
	STAR --runMode genomeGenerate --genomeDir Reference/DmelDmoj \
		--genomeFastaFiles Reference/AAA/melmoj.fa \
		--sjdbGTFfile Reference/AAA/melmoj.gtf

Reference/DmelScer:
	mkdir $@

Reference/DmelDper:
	bowtie2-build Reference/AAA/melper.fa $@
	mkdir $@

Reference/DmelDwil:
	bowtie2-build Reference/AAA/melwil.fa $@
	mkdir $@
Reference/DmelDvir:
	bowtie2-build Reference/AAA/melvir.fa $@
	mkdir $@
Reference/DmelDmoj:
	bowtie2-build Reference/AAA/melmoj.fa $@
	mkdir $@
