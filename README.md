# Germline variant discovery
## Contents
- [Data preprocessing](#Data-preprocessing)
    - [Quality control](#Quality-control)
    - [Trimming](#Trimming)
    - [Alignment](#Alignment)
    - [Lane merging (optional)](#Lane-merging-(optional))
    - [Cleaning up alignments](#Cleaning-up-alignments)
       - [Indel realignment (deprecated now)](#Indel-realignment-(deprecated-now))
       - [Mark duplicates](#Mark-duplicates)
       - [Base Quality Scores Recalibration](#Base-Quality-Scores-Recalibration)
    - [Extract metrics](#Extract-metrics)
- [Variant discovery](#Variant-discovery)
- [Callset refinement](#Callset-refinement)
   - [Hard filtering](#Hard-filtering)
   - [VQSR](#VQSR)
   - [CNN](#CNN)
- [Evaluate callset](#Evaluate-callset)
- [Functional Annotation](#Functional-Annotation)
 
![](https://github.com/R-K94/Germline_variant_discovery/blob/main/01.png)



## Data preprocessing <a name="Data-preprocessing"></a>
### Purpose
The is the obligatory first phase that must precede all variant discovery. It involves pre-processing the raw sequence data (provided in FASTQ format) to produce analysis-ready BAM files. This involves alignment to a reference genome as well as some data cleanup operations to correct for technical biases and make the data suitable for analysis.

### Input
FASTQ format, Compressed or uncompressed, single end or paired end.

### Output
BAM file and its index file ready for variant discovery

### Tools
[BWA](https://github.com/lh3/bwa), [SAMtools](https://www.htslib.org/), [Bcftools](https://www.htslib.org/download/), [GATK](https://github.com/broadinstitute/gatk), [multiqc](https://multiqc.info/)
A way to save space while working is to pipe the commands together. **Pipe tutorial [here](https://www.biostars.org/p/43677/)**.

### Required data
GATK needs a bunch of databases and has constructed several series of reference files stored in [GATK bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

Download the following from GATK bundle:

#### Truth data sets
```
1000G_phase1.snps.high_confidence.hg38.vcf
1000G_omni2.5.hg38.vcf
hapmap_3.3.hg38.vcf
Mills_and_1000G_gold_standard.indels.hg38.vcf
af-only-gnomad.hg38.vcf.gz
HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
resources_broad_hg38_v0_Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf
```

#### Reference genome
The latest version of hg38 reference genome can be downloaded from [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/),[ENSEMBL](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/) or [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/).

There are some descrepancies between these genomes and NCBI version is recommended ([see the details here.](https://genome.ucsc.edu/FAQ/FAQgenes.html#ncbiRefseq)).
Also, NCBI uses different chromosome naming than UCSC and ENSEMBL. The namings can be seen in ```*_assembly_report.txt``` file in the genome directory. For example the information on GRCh38.p13 is [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt).
- The NCBI nomenclature is like this: NC_000001.11
- The ensembl nomenclature is like this: 1
- The UCSC nomenclature is like this: chr


#### dbSNP
The latest version of dbsnp can be downloaded from [here](https://ftp.ncbi.nih.gov/snp/latest_release/).
GATk has its version of dbsnp in their bundle but its old.

### Notes
- These databases use diffrent chromosome nomenclatures which could cause problems. So, we need to make sure all the databases use the same nomenclature. 
#You can convert chr names with ```bcftools```
```
#if it is zipped first unzip:
gzip -d file.vcf.gz

#then
bcftools annotate --rename-chrs namesi.tsv yourfile.vcf > yourfile2.vcf
```
* ```namesi.tsv``` is a table that has the old and new names in two columns.

- To zip the files use the following command:
```
bgzip -c yourfile.vcf > yourfile.vcf.gz

#index
tabix -f -p vcf gzip -d file.vcf.gz
```


-----------------------------------------------------
### Steps
#### 1. Quality control <a name="Quality-control"></a>

Check the files
```
cat R1.fastq | head -n4
cat R2.fastq | head -n4
```

You could also count the reads:
```
zgrep -c "^@+" R1.fastq.gz
```

Check the quality
```
Fastqc yourfiles
```

#### 2. Trimming <a name="Trimming"></a>

At this stage, we need to trim bad reads particularly **bad 3’ ends** and **adapter sequences**.

Although nowadays this doesn’t happen often, it does still happen. In some cases, miRNA, it is expected to have adapters. Since they are not part of the genome of interest they should be removed if enough reads have them.
Also, Poly-G artifacts appear on two-channel sequencing system when the dark base G is called after synthesis has terminated. This results in calling several erroneous high-confidence G bases on the ends of affected reads. For contaminated samples, a large number of affected reads can be mapped to reference regions with high G content. This can cause problems for processing downstream.
 
To be able to remove adapters and low qualtity bases, we will use Trimmomatic.
The adapter files is intrimmomatic folder.

```
java -jar trimmomatic-0.39.jar PE \
input_forward.fq.gz input_reverse.fq.gz \
output_forward_paired.fq.gz \
output_forward_unpaired.fq.gz \
output_reverse_paired.fq.gz \
output_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads TRAILING:20 MINLEN:32 \
-trimlog /Results/log.txt
```

```-trimlog /Results/log.txt``` makes a log file to be used in ```multiqc``` and may result in a huge file. Omit this part if you don't have enough space.

***Trimmomatic Options:***
- PE: paired end
- ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
- SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
- LEADING: Cut bases off the start of a read, if below a threshold quality
- TRAILING: Cut bases off the end of a read, if below a threshold quality
- CROP: Cut the read to a specified length
- HEADCROP: Cut the specified number of bases from the start of the read
- MINLEN: Drop the read if it is below a specified length
- TOPHRED33: Convert quality scores to Phred-33
- TOPHRED64: Convert quality scores to Phred-64

Check the quality again
```
Fastqc output_forward_paired files
```

***Optional:*** Aggregate all quality and trimming data
```
Multiqc yourfilesdirectory
```

#### 3. Alignment <a name="Alignment"></a>


**Index reference** 
Before alignment Reference genome index and dictionary is needed for alingment and haplotype calling.
```
bwa index -p prefix reference.fa

gatk CreateSequenceDictionary -R reference.fa
samtools faidx reference.fa
```

After trimming, the raw reads are cleaned up of artefacts so we can align the read to the reference.
```
bwa mem -M -t 4 \
      
  -R '@RG\tID:'${sample}'\tSM:'${sample}'\tLB:'${sample}'\tPU:run${sample}_1\tCN:Your Institute\tPL:ILLUMINA' \
  ${Reference_genome} \
  paired_1.fq.gz \
  paired_2.fq.gz \
  | gatk SortSam \
  -I /dev/stdin \
  -O ${sample}.sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000
  ```
  
 ``` -R '@RG\tID:'${sample}'\tSM:'${sample}'\tLB:'${sample}'\tPL:Illumina'``` means adding read group label to output file, which is required by GATK. ```AddOrReplaceReadGroups``` in GATK can the same thing, either. The read group in BAM file looks like:```@RG	ID:test	SM:test	LB:test	PL:Illumina```



***optional:*** You can check bam structure:
 
```samtools view ${sample}.sorted.bam| head -n4```

![](https://us.v-cdn.net/5019796/uploads/editor/f4/uuzmf2cbau1y.png)

- You can get additional info on sam files: 
[samtools github](https://samtools.github.io/hts-specs/SAMtags.pdf), [samtools manual](https://samtools.sourceforge.net/SAM1.pdf), [quinlanlab](http://quinlanlab.org/tutorials/samtools/samtools.html#:~:text=Indexing%20a%20genome%20sorted%20BAM,region%20to%20which%20you%20navigate.&text=This%20will%20create%20an%20additional%20%E2%80%9Cindex%E2%80%9D%20file).
- An in depth explanation of the CIGAR can be found [here](https://genome.sph.umich.edu/wiki/SAM). The exact details of the cigar string can be found in the SAM spec as well.



***Some Additional notes***
- Say you want to count the *un-aligned* reads, you can use
```samtools view -c -f4 ${sample}.sorted.bam```
- Or you want to count the *aligned* reads you, can use
```samtools view -c -F4 ${sample}.sorted.bam```
- To subset from a bam file, first make bam and sort it then subset then sort again.
```
samtools view -h ${sample}.sorted.bam NC_000001.11:17000000-18678365 | \
gatk SortSam \
-I /dev/stdin \
  -O ./chr1.sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000
  ```
  
- if you want to convert bam to fastq:  
```bedtools bamtofastq -i ./chr1.sorted.bam \
                      -fq ./chr1.end1.fq \
                      -fq2 ./chr1.end2.fq```


#### 4. Lane merging (optional) <a name="Lane-merging-(optional)"></a>
#### 5. Cleaning up alignments <a name="Cleaning-up-alignments"></a>
 - 5.1 
 - 5.2
 - 5.3
#### 6. Extract metrics <a name="Trimming"></a>

## Variant discovery  <a name="Variant-discovery"></a>
### Purpose
### Input
### Output
### Tools
### References
---------------------------------------------------------
### Steps

## Callset refinement <a name="Callset-refinement"></a>
### Purpose
### Input
### Output
### Tools
### References
------------------------------------------------------------
### Steps

## Evaluate callset <a name="Evaluate-callset"></a>

## Functional annotation <a name="Functional-annotation"></a>

