# Germline variant discovery
## Contents
- [Data preprocessing](#Data-preprocessing)
    - [Quality control](#Quality-control)
    - [Trimming](#Trimming)
    - [Alignment](#Alignment)
    - [Lane merging (for multiple samples)](#Lane-merging-(optional))
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
GATK needs a bunch of databases and has constructed several series of reference files stored in [GATK bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) and [gatk best practices](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false).

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

## Notes<a name="Notes"></a>
#### converting database nomenclature
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
#### Intervals
For exome sequencing datasets you should always use an interval file. As explained on [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852?id=11009-):
>For exomes and similarly targeted data types, the interval list should correspond to the capture targets used for the library prep, and is typically provided by the prep kit manufacturer (with versions for each ref genome build of course). When we use intervals in our production pipelines with targeted sequencing, we make sure to give sufficient padding around the targeted sites (100 bp on each side). 

***In a nutshell***
- **Whole genome analysis:**
Intervals are not required but they can help speed up analysis by eliminating "difficult" regions and enabling parallelism
- **Exome analysis and other targeted sequencing:**
You must provide the list of targets, with padding, to exclude off-target noise. This will also speed up analysis and enable [parallelism](https://gatk.broadinstitute.org/hc/en-us/articles/360035532012).


Thus you should at each step where an ```-L``` option is available use the corresponding interval file. The exome kit vendor should provide you with the corresponding interval file. If not you can use gencode :

```
#Method 1

#Download latest gencode annotation from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
zcat gencode.v31.annotation.gtf.gz | awk '$3 == "exon" { print $1, $4, $5, $7, $18}' | tr ' ' '\t' | tr -d '";' >> interval_list.bed

#Method 2
$ curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz" |\
    gunzip -c | cut -f 3,5,6 | sort -t $'\t' -k1,1 -k2,2n | bedtools merge -i - > exome.bed
    
    
#Remove empty contigs in R
a <- read.delim("/media/myassi/01CD77B76BF0B4F0/NGS/RK/databases/Gencode/interval_list.bed", header = F)

a$V6 <- ifelse(a$V2 < a$V3, "YES", "NO")
a[rownames(a[a$V6 == "NO",] ),]

d <- a[rownames(a[a$V6 == "YES",] ),-6]

b <-  read.delim("/media/myassi/01CD77B76BF0B4F0/NGS/RK/databases/Gencode/ucsctoncbi.csv", header = F,
                 sep = ",")
b$N <- seq(1,length(b$V1))
c <- merge(d,b, all = TRUE, by.x= "V1", by.y = "V2")
e <- c[,c(6,2,3,4,5)]
write.table(e,"/media/myassi/01CD77B76BF0B4F0/NGS/RK/databases/Gencode/interval_list.bed", quote = F, row.names = F, col.names = F, sep = "\t")


# convert to interval list file
gatk BedToIntervalList -I interval_list.bed -O \ 
list.interval_list -SD GRCh38.dict

```


For information if you use GATK4 best practices these are the tools where the interval file needs to be specified
- BaseRecalibrator
- ApplyBQSR
- HaplotypeCaller
- GenomicsDBimport
- GenotypeGVCFs
- VariantRecalibrator
- ApplyVQSR


some useful links for more information: [1](https://www.biostars.org/p/486600/), [2](https://www.biostars.org/p/5187/), [3](https://www.biostars.org/p/273234/), [4](https://www.biostars.org/p/486600/)

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
```
bedtools bamtofastq -i ./chr1.sorted.bam \
                      -fq ./chr1.end1.fq \
                      -fq2 ./chr1.end2.fq
```


#### 4. Lane merging (for multiple samples) <a name="Lane-merging-(optional)"></a>
In case we generate multiple lane of sequencing or mutliple library. It is not practical to keep the data splited and all the reads should be merge into one massive file.
Since we identified the reads in the BAM with read groups, even after the merging, we can still identify the origin of each read.

```
gatk MergeSamFiles \
-AS false \
—-CREATE_INDEX true \                 
-I <input.bam>]  \
-MSD false \
-O <output_path> \
-SO coordinate \
—-USE_THREADING true \
—-VALIDATION_STRINGENCY STRICT
```

#### 5. Cleaning up alignments <a name="Cleaning-up-alignments"></a>

**5.1 Indel realignment (deprecated now)**<a name="Indel-realignment-(deprecated-now)"></a>
 
 The first step for this is to realign around indels and snp dense regions. The Genome Analysis toolkit has a tool for this called IndelRealigner.
It basically runs in 2 steps:
1. Find the targets
2. Realign them

for more information visit [this link](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md).
 
**5.2 Mark duplicates**<a name="Mark-duplicates"></a>

This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates. [More info](https://gatk.broadinstitute.org/hc/en-us/articles/360036359852-MarkDuplicates-Picard-)

```
gatk MarkDuplicates \
     -I input.bam \
     -O marked_duplicates.bam \
     -M marked_dup.metrics
```

**5.3 Base Quality Scores Recalibration** <a name="Base-Quality-Scores-Recalibration"></a>

Before performing this step, read the [Notes](#Notes) carefully.
The goal for this step is to try to recalibrate base quality scores. The vendors tend to inflate the values of the bases in the reads. Also, this step tries to lower the scores of some biased motifs for some technologies.

It runs in 2 steps:
1. Build covariates based on context and known snp sites
2. Correct the reads based on these metrics

But before that you need to make gatk index:

Most GATK tools additionally require that the main FASTA file be accompanied by a dictionary file ending in ```.dict``` and an index file ending in ```.fai```, because it allows efficient random access to the reference bases. GATK will look for these index files based on their name, so it is important that they have the same basename as the FASTA file. If you do not have these files available for your organism's reference file, you can generate them very easily; instructions are included below.

```
gatk CreateSequenceDictionary -R reference.fa
samtools faidx reference.fa
```

##### step 1: BaseRecalibrator
```
gatk BaseRecalibrator \
   -I ${bam}.sorted.dup.bam \
   -R ${genome} \
   --known-sites ${dbsnp} \
   --known-sites ${mills} \
   --known-sites ${known_indels} \
   -O ${bam}.sorted.dup.recalibration_report.table \
   -L ${intervals} \
   -ip 100 \
```

##### step 2: ApplyBQSR
```
gatk ApplyBQSR \
  -R ${genome} \
  -bqsr ${bam}.sorted.dup.recalibration_report.table \
  -I ${bam}.sorted.dup.bam  \
  -O ${bam}.sorted.dup.recal.bam 
   -L ${intervals} \
   -ip 100 \
   ```

#### 6. Extract metrics <a name="Extract-metrics"></a>

Once your whole bam is generated, it’s always a good thing to check the data again to see if everything makes sense.

***Coverage***

```
gatk DepthOfCoverage \
  --omit-depth-output-at-each-base \
  --summary-coverage-threshold 10 \
  --summary-coverage-threshold 25 \
  --summary-coverage-threshold 50 \
  --summary-coverage-threshold 100 \
  --start 1 --stop 500 --nBins 499  \
  -R ${genome} \
  -O ${bam}.sorted.dup.recal.coverage \
  -I ${bam}.sorted.dup.recal.bam
  -L ${intervals} \
   -ip 100 \
```
```summaryCoverageThreshold``` is a usefull function to see if your coverage is uniform. Another way is to compare the mean to the median. If both are almost equal, your coverage is pretty flat. If both are quite different, that means something is wrong in your coverage.


Look at the coverage

```
less -S ${bam}.sorted.dup.recal.coverage.sample_interval_summary
```

***Insert size***

```
gatk CollectInsertSizeMetrics \
  -R ${genome} \
  -I ${bam}.sorted.dup.recal.bam \
  -O ${bam}.sorted.dup.recal.metric.insertSize.tsv \
  -H ${bam}.sorted.dup.recal.metric.insertSize.histo.pdf \
  --METRIC_ACCUMULATION_LEVEL LIBRARY

head -9 ${bam}.sorted.dup.recal.metric.insertSize.tsv | tail -3 | cut -f6,7

```
For more information visit [1](https://gatk.broadinstitute.org/hc/en-us/articles/360041851491-DepthOfCoverage-BETA-) and [2](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-)

The final size of the fragment is not a big issue if it's generally longer than the sequencing design (here 200bp) to avoid overlapp. That being said, longer fragment basically could help to catch larger SV events.
The precision is important because smaller SD will allow detecting a wider range of SV events.

***Alignment metrics***

```
CollectAlignmentSummaryMetrics \
         R=reference_sequence.fasta \
         I=input.bam \
         O=output.txt
--METRIC_ACCUMULATION_LEVEL LIBRARY
```

What is the percent of aligned reads ?

```
head -10  ${bam}.sorted.dup.recal.metric.alignment.tsv | tail -4 | cut -f7
```

Usually, we consider:
- A good alignment if > 90%
- Reference assembly issues if [75-90]%
- Probably a mismatch between sample and reference if < 75 %

## Variant discovery  <a name="Variant-discovery"></a>
### Purpose
The main variant caller in GATK is HaplotypeCaller which has two modes, one for single sample and another for cohort sample. It's quite easy to select mode, if your have only one sample, use single sample mode, if not, use cohort mode (Joint call). Above is the new pipeline from GATK developed for single sample which involves deep learning is variants qaulity control.

### Input

Analysis-Ready Reads (BAM format as well as its index, output of pre-processing)

### Output

A variant information file (VCF) contains SNPs and Indels, along with its index

### Tools

GATK

---------------------------------------------------------
### Steps

```
gatk HaplotypeCaller \
-R ${genome} \
-I ${bam}.sorted.dup.recal.bam \
-O ${vcf}.vcf \
-bamout ${bam}.relaligned.bam \
--tmp-dir ${TMP}
```

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

