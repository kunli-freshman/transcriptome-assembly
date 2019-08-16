# Transcriptome-assembly
Transcriptome assembly of Japanese eels (Anguilla japonica) <br>
# Original data
four files in fastq format which have been removed adapters: `FW_1.fq`,`FW_2.fq`,`SW_1.fq`,`SW_2.fq`.<br>
>*FW: freshwater, SW: seawater
# Step1: Quality control
```
fastqc -o fastqc_results *.fastq
```
#### *Required parameters
*.fastq: The sequence files you want to get quality control reports
#### *Optional parameters
-o --outdir: Create all output files in the specified output directory.<br>
-t --threads: Specifies the number of files which can be processed simultaneously.<br>
-q --quiet: Supress all progress messages on stdout and only report errors.
```
multiqc fastqc_results/ -o multiqc_results
```
#### *Required parameters
fastqc_results/: The diretory that you store your fastqc results.
#### *Optional parameters
-o, --outdir: Create report in the specified output directory.<br>
-q, --quiet: Only show log warnings.
# Step2: Transcriptome assembly
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
### Function2: 
