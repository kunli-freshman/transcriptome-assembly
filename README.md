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
TransDecoder.Predict -t target_hisat2_cufflinks
_FW.fasta --retain_pfam_hits hisat2_cufflinks
_FW_pfam.domtblout --retain_blastp_hits hisat2_cufflinks
_FW_blastp.outfmt6
```
#### *Required parameters
-t: the target file in fasta format
#### *Optional parameters
--retain_long_orfs_mode: 'dynamic' or 'strict', default: dynamic
--retain_pfam_hits <string>: domain table output file from running hmmscan to search Pfam
--retain_blastp_hits: blastp output in '-outfmt 6' format.
## 3.2 Merge CDS
```
cat target_Trinity_FW.fasta target_spades_FW.fasta arget_hisat2_trinity_FW.fasta target_hisat2_cufflinks_FW.fasta > FW.fasta
```
## 3.3 Cluster CDS: cd-hit
```
cd-hit-est -i FW.fasta -o eel -c 0.99 -n 10 -M 16000 - T 8
```
#### *Required parameters
-i: input filename in fasta format, can be in .gz format
-o: output filename
#### *Optional parameters
 -c: sequence identity threshold, default 0.9
 -M: memory limit (in MB) for the program, default 800
 -T: number of threads, default 1
