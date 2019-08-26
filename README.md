# Transcriptome-assembly
Transcriptome assembly of Japanese eels (Anguilla japonica) <br>
# Original data
four files in fastq format which have been removed adapters: `FW_1.fq`,`FW_2.fq`,`SW_1.fq`,`SW_2.fq`.<br>
>*FW: freshwater, SW: seawater
# Step1: Quality control
## 1.1 fastQC
```
fastqc -o fastqc_results *.fastq
```
#### *Required parameters
*.fastq: The sequence files you want to get quality control reports
#### *Optional parameters
-o --outdir: Create all output files in the specified output directory.<br>
-t --threads: Specifies the number of files which can be processed simultaneously.<br>
-q --quiet: Supress all progress messages on stdout and only report errors.
## 1.2 multiQC
```
multiqc fastqc_results/ -o multiqc_results
```
#### *Required parameters
fastqc_results/: The diretory that you store your fastqc results.
#### *Optional parameters
-o, --outdir: Create report in the specified output directory.<br>
-q, --quiet: Only show log warnings.
# Step2: Transcriptome assembly (show FW as exsample)
## Method1: de novo assembly
### Function1: Trinity
```
Trinity --seqType fq --max_memory 200G --CPU 8 --left FW_1.fq --right FW_2.fq --output trinity
```
#### *Required parameters
--seqType: type of reads: ('fa' or 'fq')<br>
--max_memory: A basic recommendation is to have ~1G of RAM per ~1M pairs of Illumina reads.<br>
 If paired reads:<br>
     --left: left reads, one or more file names (separated by commas, no spaces)<br>
     --right: right reads, one or more file names (separated by commas, no spaces)<br>
 Or, if unpaired reads:<br>
     --single: single reads, one or more file names, comma-delimited<br>
 Or,<br>
      --samples_file: tab-delimited text file indicating biological replicate relationships.<br>
                                  ex.<br>
                                       cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq<br>
                                       cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq<br>
                                       cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq<br>
                                       cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq<br>
                          # if single-end instead of paired-end, then leave the 4th column above empty.<br>
#### *Optional parameters
--SS_lib_type: Strand-specific RNA-Seq read orientation.if paired: RF or FR,if single: F or R. <br>
--CPU: number of CPUs to use, default: 2<br>
--min_contig_length: minimum assembled contig length to report, default= 200 <br>
--output: name of directory for output,default: your current working directory; note: must include 'trinity' in the name as a safety precaution!<br>
### Function2: SPAdes
```
spades.py --rna -1 FW_1.fq -2 FW_2.fq -o spades
```
#### *Required parameters
-1: file with forward paired-end reads
-2: file with reverse paired-end reads
-o: directory to store all the resulting files
--rna: this flag is required for RNA-Seq data
#### *Optional parameters
-t/--threads: number of threads, default: 16
-m/--memory: RAM limit for SPAdes in Gb (terminates if exceeded), default: 250
-k: comma-separated list of k-mer sizes (must be odd and less than 128), default: 'auto'
## Method2: reference_based assembly
### Fuction1: hisat2 + trinity
### hisat2: build index
```
hisat2-build -p 4 reference/eel.fa hisat2_index/eel_tran
```
#### *Required parameters
reference_in: 'reference/eel.fa', comma-separated list of files with ref sequences
hisat2_index_base: 'hisat2_index/eel_tran', write ht2 data to files with this dir/basename
#### *Optional parameters
-p: number of threads
--exon: Exon file name
--ss: Splice site file name
### hisat2: mapping
```
hisat2 -t -x hisat2_index/eel_tran -1 FW_1.fq -2 FW_2.fq -S hisat2_FW.sam
```
#### *Required parameters
ht2-idx: 'hisat2_index/eel-tran', Index filename prefix (minus trailing .X.ht2)<br>
If paired reads:<br>
     -1: Files with #1 mates<br>
     -2: Files with #2 mates<br>
 Or, if unpaired reads:<br>
     -U: Files with unpaired reads.<br>
-S: 'hisat2_FW.sam', file for sam output, default: stdout
#### *Optional parameters
--min-intronlen: minimum intron length, default: 20<br>
--max-intronlen: maximum intron length, default: 500000<br>
 -k: the number of alns which per read report up to, default: 5<br>
 -p/--threads: number of alignment threads to launch, default: 1<br>
### samtools: convert .sam to .bam
```
samtools sort -o hisat2_FW.bam  hisat2_FW.sam
```
#### *Required parameters
hisat2_FW.sam: the file which you want convert to .bam
#### *Optional parameters
-o: the file which will be written final output rather than standard output<br>
-@, --threads: Number of additional threads to use, default: 0
### trinity
```
Trinity --genome_guided_bam hisat2_FW.bam --CPU 6 --max_memory 50G --genome_guided_max_intron 10000 --output hisat2_trinity_FW
```
#### *Required parameters
--genome_guided_bam: provide path to coordinate-sorted bam file, the file got after 'samtools: convert .sam to .bam'
#### *Optional parameters
--CPU: number of CPUs to use, default: 2<br>
--max_memory: suggested max memory to use by Trinity where limiting can be enabled<br>
--genome_guided_max_intron: use a maximum intron length that makes most sense given your targeted organism <br> 
--output: tthe directory in which the output is stored, default is a trinity_out_dir/ in your new workspace, and in this case it'll contain the resulting assembly as 'trinity-GG.fasta'. <br>
### Fuction2: hisat2 + cufflinks
### cufflinks
```
cufflinks -p 8 -o hisat2_cufflinks_FW  hisat2_FW.bam
```
#### *Required parameters
hisat2_FW.bam: the result of 'samtools: convert .sam to .bam'
#### *Optional parameters
-o/--output-dir: write all output files to this directory, default: ./   <br>
-p/--num-threads: number of threads used during analysis , default: 1     <br>
--seed: value of random number generator seed , default: 0   <br>
-g/--GTF-guide: use reference transcript annotation to guide assembly <br>
# Step3: Finding CDS (show FW as exsample)
## 3.1 Finding CDS
### 3.1.0 convert .gff to .fasta: Transdecoder
```
gtf_genome_to_cdna_fasta.pl hisat2_cufflinks_FW.gtf reference.fasta > hisat2_cufflinks_FW.fasta
```
#### *Required parameters
hisat2_cufflinks_FW.gtf: the result of 'hisat2 + cufflinks'<br>
reference.fasta: the reference of Japanese eel<br>
hisat2_cufflinks_FW.fasta: the output of this work<br>
### 3.1.1 select the signal best ORF: Transdecoder
```
TransDecoder.LongOrfs -t hisat2_cufflinks_FW.fasta
TransDecoder.LongOrfs -t hisat2_trinity_FW.fasta
TransDecoder.LongOrfs -t trinity_FW.fasta
TransDecoder.LongOrfs -t spades_FW.fasta
```
#### *Required parameters
.fasta: the file in fasta in format which you get from assembly
#### *Optional parameters
--output_dir | -O: path to intended output directory, default: basename( -t val ) + ".transdecoder_dir"<br>
### 3.1.2 identify ORFs with homology to known proteins: blastp
```
blastp -query Trinity_FW.fasta.transdecoder_dir/longest_orfs.pep -db database/uniprot_sprot -outfmt 6 -evalue 1e-5 -num_threads 10 > Trinity_FW.fasta_blastp.outfmt6
blastp -query spades_FW.fasta.transdecoder_dir/longest_orfs.pep -db database/uniprot_sprot -outfmt 6 -evalue 1e-5 -num_threads 10 > spades_FW.fasta_blastp.outfmt6
blastp -query hisat2_trinity_FW.fasta.transdecoder_dir/longest_orfs.pep -db database/uniprot_sprot -outfmt 6 -evalue 1e-5 -num_threads 10 > hisat2_trinity_FW.fasta_blastp.outfmt6
blastp -query hisat2_cufflinks_FW.fasta.transdecoder_dir/longest_orfs.pep -db database/uniprot_sprot -outfmt 6 -evalue 1e-5 -num_threads 10 > hisat2_cufflinks_FW.fasta_blastp.outfmt6
```
#### *Required parameters
-query: the file used to align to database which is the result of 'select the signal best ORF: Transdecoder'
-db: BLAST database name<br>
-outfmt: alignment view options: 6 = Tabula<br>
#### *Optional parameters
-max_target_seqs: Maximum number of aligned sequences to keep, default = 500<br>
-num_threads: Number of threads (CPUs) to use in the BLAST search, default = 1<br>
-evalue: Expectation value (E) threshold for saving hits, default = `10'
*** Consider of the time and memory, I use Diamond to replace blastp, just as
```
diamond blastp -q Trinity_FW.fasta.transdecoder_dir/longest_orfs.pep -d database/uniprot_sprot --outfmt 6 --evalue 1e-5 -o Trinity_FW.fasta_blastp.outfmt6
```
### 3.1.3 search for protein signature: hmmer
```
hmmscan --cpu 8 --domtblout Trinity_FW_pfam.domtblout database/Pfam-A.hmm Trinity_FW.fasta.transdecoder_dir/longest_orfs.pep
hmmscan --cpu 8 --domtblout spades_FW_pfam.domtblout database/Pfam-A.hmm spades_FW.fasta.transdecoder_dir/longest_orfs.pep
hmmscan --cpu 8 --domtblout hisat2_trinity_FW_pfam.domtblout database/Pfam-A.hmm hisat2_trinity_FW.fasta.transdecoder_dir/longest_orfs.pep
hmmscan --cpu 8 --domtblout hisat2_cufflinks_FW_pfam.domtblout database/Pfam-A.hmm hisat2_cufflinks_FW.fasta.transdecoder_dir/longest_orfs.pep
```
#### *Required parameters
hmmdb: Pfam database name,'database/Pfam-A.hmm'<br>
seqfile: the file used to align to database which is the result of 'select the signal best ORF: Transdecoder'
#### *Optional parameters
--domtblout: save parseable table of per-domain hits to file<br>
--cpu: number of parallel CPU workers to use for multithreads, default = 2
### 3.1.4 predict CDS: Transdecoder
```
TransDecoder.Predict -t target_Trinity_FW.fasta --retain_pfam_hits Trinity_FW_pfam.domtblout --retain_blastp_hits Trinity_FW_blastp.outfmt6
TransDecoder.Predict -t target_spades_FW.fasta --retain_pfam_hits spades_FW_pfam.domtblout --retain_blastp_hits spades_FW_blastp.outfmt6
TransDecoder.Predict -t target_hisat2_trinity_FW.fasta --retain_pfam_hits hisat2_trinity_FW_pfam.domtblout --retain_blastp_hits hisat2_trinity_FW_blastp.outfmt6
TransDecoder.Predict -t target_hisat2_cufflinks_FW.fasta --retain_pfam_hits hisat2_cufflinks_FW_pfam.domtblout --retain_blastp_hits hisat2_cufflinks_FW_blastp.outfmt6
```
#### *Required parameters
 -t: the target file in fasta format
#### *Optional parameters
 --retain_long_orfs_mode: 'dynamic' or 'strict', default: dynamic<br>
 --retain_pfam_hits <string>: domain table output file from running hmmscan to search Pfam<br>
 --retain_blastp_hits: blastp output in '-outfmt 6' format.
## 3.2 Merge CDS
```
cat target_Trinity_FW.fasta target_spades_FW.fasta arget_hisat2_trinity_FW.fasta target_hisat2_cufflinks_FW.fasta > FW.fasta
```
## 3.3 Cluster CDS: cd-hit
```
cd-hit-est -i FW.fasta -o eel -c 0.99 -M 16000 -T 8
```
#### *Required parameters
 -i: input filename in fasta format, can be in .gz format<br>
 -o: output filename
#### *Optional parameters
 -c: sequence identity threshold, default 0.9<br>
 -M: memory limit (in MB) for the program, default 800<br>
 -T: number of threads, default 1<br>
 -n：word_length, default 10, see user's guide for choosing it<br>
# Step4: Transcriptome annotation & gene ontology (show FW as exsample)
## 4.1 Transcriptome annotation: blast
### 4.1.0 make taxidlist of NR, NT, Swissprot(Japanese eel belongs to Metazoa)
```
 get_species_taxids.sh  -n Metazoa
 get_species_taxids.sh  -t 33208 > list
 ```
### 4.1.1 compare to NR
```
 blastx -query cd_hit/eel_FW -db database/nr -taxidlist list -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out FW.nr.outfmt6
```
### 4.1.2 compare to NT
 ```
 blastn -query cd_hit/eel_FW -db database/nt -taxidlist list -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out FW.nt.outfmt6
```
### 4.1.3 compare to Swissprot
```
blastx -query cd_hit/eel_FW -db database/swissprot -taxidlist list -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out FW.swissprot.outfmt6
```
#### *Required parameters
-query: Input file name<br>
-db: BLAST database name<br>
#### *Optional parameters
-out: Output file name, default = `-'<br>
-outfmt: alignment view options, 5 = BLAST XML<br>
-evalue: Expectation value (E) threshold for saving hits, default = 10
-num_threads: Number of threads (CPUs) to use in the BLAST search, default = 1
## 4.2 Gene ontology: blast2GO
```
???? java -jar blast2go.jar -in FW.nr.outfmt5 -a -v -out FW_GO_nr
```
# Step5: Transcriptome assembly validation (show FW as exsample)
## 5.1 quality assessment: rnaQUAST
 ```
python rnaQUAST.py --threads 8 --transcripts cd_hit/eel_FW -1 FW_1.fq -2 FW_2.fq --busco_lineage busco/actinopterygii_odb9 -o rnaQUAST/FW_rnaQUAST
 ```
#### *Required parameters
-1, --left_reads: File with forward paired-end reads in FASTQ or gzip-compressed fastq format.<br>
-2, --right_reads: File with reverse paired-end reads in FASTQ or gzip-compressed fastq format. <br>
-o, --output_dir: Directory to store all results. Default is rnaQUAST_results/results_<datetime>. <br>
#### *Optional parameters
--prokaryote: Use this option if the genome is prokaryotic. <br>
-ss, --strand_specific: Set if transcripts were assembled using strand-specific RNA-Seq data in order to benefit from knowing whether the transcript originated from the + or - strand. <br>
--min_alignment: Minimal alignment length to be used, default value is 50. <br>
--busco_lineage: Run BUSCO tool, which detects core genes in the assembly.<br>
--tophat: Run with TopHat tool instead of STAR for analyzing database coverage by reads.<br>
-t, --threads: Maximum number of threads. Default is min(number of CPUs / 2, 16). <br>
## 5.2 quality control: detonate
 ```
 rsem-eval-estimate-transcript-length-distribution cd-hit/eel_FW denotate/parameter_file
 ```
 ```
rsem-eval-calculate-score --transcript-length-parameters denotate/parameter_file --paired-end FW_1.fq FW_2.fq cd_hit/eel_FW eel_FW 737
 ```
#### *Required parameters
--paired-end upstream_read_file(s) downstream_read_file(s): "FW_1.fq FW_2.fq", the original data in fastq format.<br>
assembly_fasta_file: "cd_hit/eel_FW", a multi-FASTA file contains the assembly used for calculating RSEM-EVAL score.<br>
sample_name: "eel_FW",  The name of the sample analyzed. All output files are prefixed by this name.<br>
L: "737",  For paired-end data, L represents the average fragment length. <br>
#### *Optional parameters
--overlap-size: The minimum overlap size required to join two reads together. Default = 0<br>
--transcript-length-parameters: Read the true transcript length distribution's mean and standard deviation from <file><br>
# Step6: difference analysis (show FW as exsample)
## 6.1 gene expression level: RSEM
```
rsem-prepare-reference --bowtie2 cd_hit/eel_FW reference/eel_FW
```
#### *Required parameters
reference_fasta_file(s): "", either a comma-separated list of Multi-FASTA formatted files OR a directory name. If a directory name is specified, RSEM will read all files with suffix ".fa" or ".fasta" in this directory. The files should contain either the sequences of transcripts or an entire genome, depending on whether the '--gtf' option is used.<br>
reference name: "", The name of the reference used. RSEM will generate several reference-related files that are prefixed by this name. This name can contain path information.
#### *Optional parameters
 --gtf: If this option is on, RSEM assumes that 'reference_fasta_file(s)' contains the sequence of a genome, and will extract transcript
        reference sequences using the gene annotations specified in <reference_fasta_file(s)>, which should be in GTF format.
        If this and '--gff3' options are off, RSEM will assume 'reference_fasta_file(s)' contains the reference transcripts. In
        this case, RSEM assumes that name of each sequence in the Multi-FASTA files is its transcript_id.<br>
 --transcript-to-gene-map <file>: Use information from <file> to map from transcript (isoform) ids to gene ids. Each line of <file> should be of the form:<br>
        gene_id transcript_id<br>
with the two fields separated by a tab character.<br>
 --bowtie: Build Bowtie indices. <br>
 --bowtie-path <path>: The path to the Bowtie executables. Default: the path to Bowtie executables is assumed to be in the user's PATH environment variable<br>
--bowtie2: Build Bowtie 2 indices. Default: off<br>
--bowtie2-path: The path to the Bowtie 2 executables. Default: the path to Bowtie 2 executables is assumed to be in the user's PATH environment variable<br>
--star: Build STAR indices. Default: off<br>
--star-path <path>: The path to STAR's executable. Default: the path to STAR executable is assumed to be in user's PATH environment variable<br><br>
  -p/--num-threads <int>: Number of threads to use for building STAR's genome indices.<br>
 ```
rsem-calculate-expression -p 8 --paired-end --bowtie2 --estimate-rspd --append-names FW_1.fq FW_2.fq reference/eel_FW exp/eel_FW
 ```
#### *Required parameters
upstream_read_files(s): Comma-separated list of files containing single-end reads or upstream reads for paired-end data. By default, these files are assumed to be in FASTQ format. If the --no-qualities option is specified, then FASTA format is expected.<br>
downstream_read_file(s): Comma-separated list of files containing downstream reads which are paired with the upstream reads. By default, these files are assumed to be in FASTQ format. If the --no-qualities option is specified, then FASTA format is expected.<br>
reference_name: The name of the reference used. The user must have run 'rsem-prepare-reference' with this reference_name before running
this program.<br>
sample_name: The name of the sample analyzed. All output files are prefixed by this name <br>
#### *Optional parameters
--paired-end: Input reads are paired-end reads. Default: off<br>
-p/--num-threads <int>: Number of threads to use. Both Bowtie/Bowtie2, expression estimation and 'samtools sort' will use this many threads. Default: 1<br>
--bowtie2: Use Bowtie 2 instead of Bowtie to align reads. Since currently RSEM does not handle indel, local and discordant alignments, the Bowtie2 parameters are set in a way to avoid those alignments. In particular, we use options '--sensitive --dpad 0 --gbar 99999999
--mp 1,1 --np 1 --score-min L,0,-0.1' by default. The last parameter of '--score-min', '-0.1', is the negative of maximum mismatch rate.
This rate can be set by option '--bowtie2-mismatch-rate'. If reads are paired-end, we additionally use options '--no-mixed' and
'--no-discordant'. Default: off<br>
--star: Use STAR to align reads. Alignment parameters are from ENCODE3's STAR-RSEM pipeline. To save computational time and memory resources, STAR's Output BAM file is unsorted. It is stored in RSEM's temporary directory with name as 'sample_name.bam'. Each STAR job will have its own private copy of the genome in memory. Default: off<br>
 --append-names: If gene_name/transcript_name is available, append it to the end of gene_id/transcript_id (separated by '_') in files
'sample_name.isoforms.results' and 'sample_name.genes.results'. Default: off<br>
 --estimate-rspd: Set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD. Default: off<br> 
## 6.2 difference analysis: EBSeq
```
rsem-generate-ngvector reference/eel eel_ref
rsem-generate-data-matrix FW.isoforms.results SW.isoforms.results > IsoMat.txt
rsem-run-ebseq --ngvector eel_ref.ngvec IsoMat.txt 1,1 IsoMat.results
rsem-control-fdr IsoMat.results 0.05 IsoMat.de.txt
```
## 6.3 GO enrichment analysis
