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

#

### Input
FASTQ format, Compressed or uncompressed, single end or paired end.

#

### Output
BAM file and its index file ready for variant discovery

#

### Tools
[BWA](https://github.com/lh3/bwa), [SAMtools](https://www.htslib.org/), [Bcftools](https://www.htslib.org/download/), [GATK](https://github.com/broadinstitute/gatk), [multiqc](https://multiqc.info/)
A way to save space while working is to pipe the commands together. **Pipe tutorial [here](https://www.biostars.org/p/43677/)**.

#

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
GATK recommends using hg38 [more details](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19).

> We strongly recommend switching to GRCh38/hg38 if you are working with human sequence data. In addition to adding many alternate contigs, GRCh38 corrects thousands of small sequencing artifacts that cause false SNPs and indels to be called when using the GRCh37 assembly (b37/Hg19). It also includes synthetic centromeric sequence and updates non-nuclear genomic sequence.

The latest version of hg38 reference genome can be downloaded from [NCBI (GCF/Refseq)](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/),[NCBI (GCA/Genebank)](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/), [ENSEMBL](https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/) or [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/).

[This link](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) also may provide useful information.

There are some descrepancies between these genomes and NCBI version is recommended ([see the details here.](https://genome.ucsc.edu/FAQ/FAQgenes.html#ncbiRefseq)).
Also, NCBI uses different chromosome naming than UCSC and ENSEMBL. The namings can be seen in ```*_assembly_report.txt``` file in the genome directory. For example the information on GRCh38.p13 is [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt).
- The NCBI nomenclature is like this: NC_000001.11 for GCF and CM000663.2 for GCA
- The ensembl nomenclature is like this: 1
- The UCSC nomenclature is like this: chr

But if you use the `GRCh38_major_release_seqs_for_alignment_pipelines` from ncbi, it has UCSC nomenclature

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
#

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

#
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


#
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

#
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

##### 5.3.1 BaseRecalibrator
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

##### 5.3.2 ApplyBQSR
```
gatk ApplyBQSR \
  -R ${genome} \
  -bqsr ${bam}.sorted.dup.recalibration_report.table \
  -I ${bam}.sorted.dup.bam  \
  -O ${bam}.sorted.dup.recal.bam 
   -L ${intervals} \
   -ip 100 \
   ```
#
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
  -I ${bam}.sorted.dup.recal.bam \
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
         O=output.txt \
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
#
### Input

Analysis-Ready Reads (BAM format as well as its index, output of pre-processing)

#
### Output

A variant information file ([VCF](https://www.internationalgenome.org/wiki/Analysis/vcf4.0)) contains SNPs and Indels, along with its index

#
### Tools

GATK HaplotypeCaller
for more info about HaplotypeCaller check [this link](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)

---------------------------------------------------------
### Steps

#### Variant Calling <a name="Variant-discovery"></a>
You can perform this step in vcf mode or in gvcf mode. For more information check [this link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format#:~:text=The%20key%20difference%20between%20a,a%20cohort%20in%20subsequent%20steps.)

**VCF calling workflow**

```
gatk HaplotypeCaller \
-R ${genome} \
-I ${bam}.sorted.dup.recal.bam \
-O ${vcf}.vcf \
-bamout ${bam}.relaligned.bam \
--tmp-dir ${TMP}
```
#
**GVCF calling workflow**
```
gatk HaplotypeCaller \
-R ${genome} \
-I ${bam}.sorted.dup.recal.bam \
-O ${vcf}.g.vcf.gz \
-bamout ${bam}.relaligned.bam \
--tmp-dir ${TMP}
-ERC GVCF
```
If you have used GVCS calling, now you need to joint call the genotypes:


#
##### joint genotyping on a singular sample

``` 
gatk GenotypeGVCFs \
   -R ${genome} \
   -V ${vcf}.g.vcf.gz \
   -O ${vcf}.vcf.gz
   
   ```
#
##### joint genotyping on a cohort

To gain more inforamtion on variant calling on a cohort check [this link](https://gatk.broadinstitute.org/hc/en-us/articles/360035890411?id=3893)
 
 **1. Data aggregation**
 
```
gatk GenomicsDBImport \
      -V data/gvcfs/mother.g.vcf.gz \
      -V data/gvcfs/father.g.vcf.gz \
      -V data/gvcfs/son.g.vcf.gz \
      --genomicsdb-workspace-path my_database \
      --tmp-dir=/path/to/large/tmp \
      -L 20
```

**2. Joint genotyping**
```
gatk GenotypeGVCFs \
   -R Homo_sapiens_assembly38.fasta \
   -V gendb://my_database \
   -O output.vcf.gz
```

## Notes 
Based on the FORMAT description, the genotype (GT) is the first information provided (separated by :).

Three different values are available:
- 0/0 : homozgote reference (with the REF columns value)
- 1/1 : homozgote variant (with the ALT columns value)
- 1/0 : heterozygote (with both REF and ALT columns values)

***How many SNPs were found?***
```grep -v "^#" file.vcf | wc -l```

***The number of indel can be computed using this command:***

```

grep -v "^#" file.vcf | awk '{ if(length($4) != length($5)) { print $0 } }' | wc -l

```

***Largest indels size can be computed by extending the previous command while printing the size of the indel:***

```

grep -v "^#" file.vcf | \
awk '{ if(length($4) != length($5)) { print sqrt((length($4) - length($5))^2) "\t"$0 } }' | \
sort -k1,1nr | head

```

#
## Callset refinement <a name="Callset-refinement"></a>
### Purpose
## 1. Variant Filteration
The GATK's variant calling tools are designed to be very lenient in order to achieve a high degree of sensitivity. This is good because it minimizes the chance of missing real variants, but it does mean that we need to filter the raw callset they produce in order to reduce the amount of false positives, which can be quite large.

As of 2022 GATK offers different approaches to site-level variant filtration. Site-level filtering involves using INFO field annotations in filtering:
- [Hard Filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471)  
- [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612)  
- [CNN](https://gatk.broadinstitute.org/hc/en-us/articles/360037226672-CNNScoreVariants)  

#### Hard Filtering
Hard-filtering consists of choosing specific thresholds for one or more annotations and throwing out any variants that have annotation values above or below the set thresholds. By annotations, we mean properties or statistics that describe for each variant e.g. what the sequence context is like around the variant site, how many reads covered it, how many reads covered each allele, what proportion of reads were in forward vs reverse orientation, and so on.

#### VQSR
The established way to filter the raw variant callset is to use variant quality score recalibration (VQSR), which uses machine learning to identify annotation profiles of variants that are likely to be real, and assigns a VQSLOD score to each variant that is much more reliable than the QUAL score calculated by the caller. In the first step of this two-step process, the program builds a model based on training variants, then applies that model to the data to assign a well-calibrated probability to each variant call. We can then use this variant quality score in the second step to filter the raw call set, thus producing a subset of calls with our desired level of quality, fine-tuned to balance specificity and sensitivity.

In order to achieve the best exome results one needs to use an exome SNP and/or indel callset with at least 30 samples. For users with experiments containing fewer exome samples there are several options to explore:
Add additional samples for variant calling, either by sequencing additional samples or using publicly available exome bams from the 1000 Genomes Project (this option is used by the Broad exome production pipeline). Be aware that you cannot simply add VCFs from the 1000 Genomes Project. You must either call variants from the original BAMs jointly with your own samples, or (better) use the reference model workflow to generate GVCFs from the original BAMs, and perform joint genotyping on those GVCFs along with your own samples' GVCFs with GenotypeGVCFs.
You can also try using the VQSR with the smaller variant callset, but experiment with argument settings (try adding ```--maxGaussians 4``` to your command line, for example). You should only do this if you are working with a non-model organism for which there are no available genomes or exomes that you can use to supplement your own cohort.

#### CNN
Single sample variant discovery uses HaplotypeCaller in its default single-sample mode to call variants in an analysis-ready BAM file. The VCF that HaplotypeCaller emits errs on the side of sensitivity, so some filtering is often desired. To filter variants first run the CNNScoreVariants tool. This tool annotates each variant with a score indicating the model's prediction of the quality of each variant. To apply filters based on those scores run the FIlterVariantTranches tool with SNP and INDEL sensitivity tranches appropriate for your task.

## 2. Genotype Refinement

While every study can benefit from increased data accuracy, this workflow is especially useful for analyses that are concerned with how many copies of each variant an individual has (e.g. in the case of loss of function) or with the transmission (or de novo origin) of a variant in a family.

If a “gold standard” dataset for SNPs is available, that can be used as a very powerful set of priors on the genotype likelihoods in your data. For analyses involving families, a pedigree file describing the relatedness of the trios in your study will provide another source of supplemental information. If neither of these applies to your data, the samples in the dataset itself can provide some degree of genotype refinement.

for more info check [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants)

#
### Input
raw variant information file (VCF) along with its index

#
### Output
Analysis ready variant information file (VCF) along with its index

#
### Tools
filteration:  
GATK [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants)  
GATK [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)  
GATK [VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator)  
GATK [ApplyVQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360037056912-ApplyVQSR)  
GATK [CNNScoreVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037226672-CNNScoreVariants)  
GATK [FilterVariantTranches](https://gatk.broadinstitute.org/hc/en-us/articles/360037227632-FilterVariantTranches)  

refinment:    
GATK [CalculateGenotypePosteriors](https://gatk.broadinstitute.org/hc/en-us/articles/360037226592-CalculateGenotypePosteriors)  
GATK [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)  
GATK [VariantAnnotator](https://gatk.broadinstitute.org/hc/en-us/articles/360037224652-VariantAnnotator-BETA-)  


------------------------------------------------------------
### Steps
### 1. INDEL and SNP separtion
For hard filtering and VQSR you need to devide your VCF to two separate files:

```
# Subset to SNPs-only callset
gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type SNP \
    -O snps.vcf.gz

# Subset to indels-only callset
gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type INDEL \
    -O indels.vcf.gz

```

#
### 2. Vairant Filteration
- [A] Hard Filtering  
- [B] VQSR
- [C] CNN

#### [A] Hard Filtering
For more information visit [this site](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471) and [this site](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2).

To filter on multiple expressions, provide each in separate expression. For INFO level annotations, the parameter is ```-filter```, which should be immediately followed by the corresponding ```–-filter-name``` label. Here we show basic hard-filtering thresholds.

To filter based on FORMAT level annotations, use ```--genotype-filter-expression``` and ```--genotype-filter-name```.

The GATK does not recommend use of compound filtering expressions, e.g. the logical || "OR". For such expressions, if a record is null for or missing a particular annotation in the expression, the tool negates the entire compound expression and so automatically passes the variant record even if it fails on one of the expressions.

```
#SNP
gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps_filtered.vcf.gz
    
    #INDEL
    gatk VariantFiltration \ 
    -V indels.vcf.gz \ 
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \ 
    -O indels_filtered.vcf.gz

```


This produces a VCF with the same variant records now annotated with filter status. Specifically, if a record passes all the filters, it receives a PASS label in the FILTER column. A record that fails a filter receives the filter name in the FILTER column, e.g. SOR3. If a record fails multiple filters, then each failing filter name appears in the FILTER column separated by semi-colons ;, e.g. ```MQRankSum-12.5;ReadPosRankSum-8````.

The filteration settings provided are GATK's recommandation. For manual filtering you can visualize your data:
First we need to extract the required data and make a table with [VariantsToTable](https://gatk.broadinstitute.org/hc/en-us/articles/360036896892-VariantsToTable).

```
gatk VariantsToTable \
-V /media/myassi/01CD77B76BF0B4F0/NGS/RK/Analysis/Variants/NA12878.hc.SNP.GATK.giab.vcf.gz \
-F CHROM -F POS -F QUAL \
-F BaseQRankSum -F MQRankSum -F ReadPosRankSum \
-F DP -F FS -F MQ -F QD -F SOR \
-GF GQ \
-O table.txt


#if you want to compare your data calls with GIAB:
#First we add annotations of giab sets to our data
gatk VariantAnnotator \
-V /media/myassi/01CD77B76BF0B4F0/NGS/RK/Analysis/Variants/NA12878.hc.SNP.GATK.vcf.gz \
-O /media/myassi/01CD77B76BF0B4F0/NGS/RK/Analysis/Variants/NA12878.hc.SNP.GATK.giab.vcf.gz \
--resource:giab /media/myassi/01CD77B76BF0B4F0/NGS/RK/databases/GIAB/NCBI_HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
--expression giab.callsets 


gatk VariantsToTable \
-V /media/myassi/01CD77B76BF0B4F0/NGS/RK/Analysis/Variants/NA12878.hc.SNP.GATK.giab.vcf.gz \
-F CHROM -F POS -F QUAL \
-F BaseQRankSum -F MQRankSum -F ReadPosRankSum \
-F DP -F FS -F MQ -F QD -F SOR \
-F giab.callsets \
-GF GQ \
-O /media/myassi/01CD77B76BF0B4F0/NGS/RK/Analysis/Variants/NA12878.hc.SNP.GATK.giab.txt
```


Next parts are in R

```
# plotting.R script loads ggplot and gridExtra libraries and defines functions to plot variant annotations 
library(ggplot2)
install.packages("gridExtra")
library(gridExtra)

require(ggplot2, quietly = TRUE)
require(gridExtra, quietly = TRUE)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# Function for making density plots of a single annotation
makeDensityPlot <- function(dataframe, xvar, split, xmin=min(dataframe[xvar], na.rm=TRUE), xmax=max(dataframe[xvar], na.rm=TRUE), alpha=0.5) {
  
  if(missing(split)) {
    return(ggplot(data=dataframe, aes_string(x=xvar)) + xlim(xmin,xmax) + geom_density() )
  }
  else {
    return(ggplot(data=dataframe, aes_string(x=xvar, fill=split)) + xlim(xmin,xmax) + geom_density(alpha=alpha) )
  }
}

# Function for making scatter plots of two annotations
makeScatterPlot <- function(dataframe, xvar, yvar, split, xmin=min(dataframe[xvar], na.rm=TRUE), xmax=max(dataframe[xvar], na.rm=TRUE), ymin=min(dataframe[yvar], na.rm=TRUE), ymax=max(dataframe[yvar], na.rm=TRUE), ptSize=1, alpha=0.6) {
  if(missing(split)) {
    return(ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + xlim(xmin,xmax) + ylim(ymin,ymax) + geom_point(size=ptSize, alpha=alpha) )
  }
  else {
    return(ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + aes_string(color=split) + xlim(xmin,xmax) + ylim(ymin,ymax) + geom_point(size=ptSize, alpha=alpha) )
  }
}

# Function for making scatter plots of two annotations with marginal density plots of each
makeScatterPlotWithMarginalDensity <- function(dataframe, xvar, yvar, split, xmin=min(dataframe[xvar], na.rm=TRUE), xmax=max(dataframe[xvar], na.rm=TRUE), ymin=min(dataframe[yvar], na.rm=TRUE), ymax=max(dataframe[yvar], na.rm=TRUE), ptSize=1, ptAlpha=0.6, fillAlpha=0.5) {
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
  if(missing(split)){
    scatter <- ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + geom_point(size=ptSize, alpha=ptAlpha) + xlim(xmin,xmax) + ylim(ymin,ymax) 
    plot_top <- ggplot(data=dataframe, aes_string(x=xvar)) + geom_density(alpha=fillAlpha) + theme(legend.position="none") + xlim(xmin,xmax) 
    plot_right <- ggplot(data=dataframe, aes_string(x=yvar)) + geom_density(alpha=fillAlpha) + coord_flip() + theme(legend.position="none") + xlim(ymin,ymax) 
  } 
  else{
    scatter <- ggplot(data=dataframe) + aes_string(x=xvar, y=yvar) + geom_point(size=ptSize, alpha=ptAlpha, aes_string(color=split)) + xlim(xmin,xmax) + ylim(ymin,ymax) 
    plot_top <- ggplot(data=dataframe, aes_string(x=xvar, fill=split)) + geom_density(alpha=fillAlpha) + theme(legend.position="none") + xlim(xmin,xmax) 
    plot_right <- ggplot(data=dataframe, aes_string(x=yvar, fill=split)) + geom_density(alpha=fillAlpha) + coord_flip() + theme(legend.position="none") + xlim(ymin,ymax) 
  }
  legend <- get_legend(scatter)
  scatter <- scatter + theme(legend.position="none")
  temp <- grid.arrange(plot_top, legend, scatter, plot_right, ncol=2, nrow=2, widths=c(4,1), heights=c(1,4))
  return(temp)
}



library(readr)
SNP.giab <- read_delim("/media/myassi/01CD77B76BF0B4F0/NGS/RK/Analysis/Variants/NA12878.hc.SNP.GATK.giab.txt","\t", 
                             escape_double = FALSE, col_types = cols(giab.callsets = col_character()), trim_ws = TRUE)
                             #Remove col_types = cols(giab.callsets = col_character()) if you didn't annotate GIAB
                             
#Below are three different plotting functions. Uncomment them one at a time and see how the graph changes as you alter them. To uncomment, simply remove the # at the beginning of the line, and put it in front of the line you no longer want to graph.

B = makeDensityPlot(SNP.giab, "QUAL") #base quality
B1 = makeDensityPlot(SNP.giab, "QUAL", xmax=1000)

#If you annotated GIAB you can add it to plot
B2 = makeDensityPlot(SNP.giab, "QUAL", xmax=3000, split="giab.callsets")


# Change up the parameters, e.g. add 'split="giab.callsets"', examine RankSums, FS and SOR
C = makeDensityPlot(SNP.giab, "QD") # quality normalized by depth
C = makeDensityPlot(SNP.giab, "QD", split="giab.callsets")
C = makeDensityPlot(SNP.giab, "BaseQRankSum", split="giab.callsets")
C = makeDensityPlot(SNP.giab, "MQRankSum", split="giab.callsets", xmax=1, xmin=-1)#pos= ref has a better quality, neg = variant has better quality
C = makeDensityPlot(SNP.giab, "ReadPosRankSum", split="giab.callsets") #
C = makeDensityPlot(SNP.giab, "FS", split="giab.callsets") #strand bias
C = makeDensityPlot(SNP.giab, "SOR", split="giab.callsets") #strand odds ratio


# When plotting two annotations, does the combination of the two tell us anything more than either did separately?
# 
#   Try adjusting the parameters.
# Substitute in other annotations. For example, the following recreates the plot on the front page of the tutorial worksheet.
# F = makeScatterPlotWithMarginalDensity(motherSNP.giab, "QUAL", "DP", split="set", xmax=10000, ymax=100, ptSize=0.5, ptAlpha=0.05)

E = makeScatterPlotWithMarginalDensity(SNP.giab, 
                                       "QD", "DP", 
                                       split="giab.callsets", 
                                       ymax=250, 
                                       ptSize=0.5, ptAlpha=0.2)

```


#
#### [B] VQSR
##### B.1 Calculate VQSLOD tranches
```
#calculate VQSLOD tranches for SNPs 
time gatk VariantRecalibrator \
--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
	-R  ${genome} \
	-V ${vcf}.vq.SNP.vcf \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
		--resource:omni,known=false,training=true,truth=true,prior=12 databases/1000Gomni2.5/NCBI_1000G_omni2.5.hg38.vcf.gz \
		--resource:1000G,known=false,training=true,truth=false,prior=10 databases/1000G_phase1_SNP/NCBI_1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=7 databases/dbsnp/GCF_000001405.39.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15 databases/Hapmap3.3/NCBI_hapmap_3.3.hg38.vcf \
	-O ${vcf}.vq.SNP.calibrate.recal \
	--tranches-file ${vcf}.vq.SNP.calibrate.cohort_snps.tranches \
	--rscript-file ${vcf}.output.plots.R \
	--tmp-dir ${TMP} \
	-L ${intervals} \
	-ip 100 \
    
    
    
#calculate VQSLOD tranches for INDELS
gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
    -V cohort_sitesonly.vcf.gz \
    --trust-all-polymorphic \
    -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
    -mode INDEL \
    --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=12:Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10:Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2:Homo_sapiens_assembly38.dbsnp138.vcf \
    -O cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    	--tmp-dir ${TMP} \
	-L ${intervals} \
	-ip 100 \
    
```

##### B.2 Filter on VQSLOD using ApplyVQSR

```
#SNP
gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V indel.recalibrated.vcf.gz \
    --recal-file ${snps_recalibration} \
    --tranches-file ${snps_tranches} \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode SNP \
    -O snp.recalibrated.vcf.gz \
        	--tmp-dir ${TMP} \
	-L ${intervals} \
	-ip 100 \
    
    
#INDEL
gatk --java-options "-Xmx5g -Xms5g" \
    ApplyVQSR \
    -V cohort_excesshet.vcf.gz \
    --recal-file cohort_indels.recal \
    --tranches-file cohort_indels.tranches \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -mode INDEL \
    -O indel.recalibrated.vcf.gz
            	--tmp-dir ${TMP} \
	-L ${intervals} \
	-ip 100 \
    

```

**99.9% is the recommended default VQSLOD cutoff for SNPs in human genomic analysis.**

#
#### [C] CNN
Annotate a VCF with scores from a Convolutional Neural Network (CNN). This tool streams variants and their reference context to a python program, which evaluates a pre-trained neural network on each variant. The default models were trained on single-sample VCFs. The default model should not be used on VCFs with annotations from joint call-sets. The neural network performs convolutions over the reference sequence surrounding the variant and combines those features with a multilayer perceptron on the variant annotations. 2D models convolve over aligned reads as well as the reference sequence, and variant annotations. 2D models require a SAM/BAM file as input and for the --tensor-type argument to be set to a tensor type which requires reads, as in the example below. Pre-trained 1D and 2D models are included in the distribution. It is possible to train your own models with the tools: CNNVariantWriteTensors and CNNVariantTrain. CNNVariantTrain will create a json architecture file and an hd5 weights file, which you can use with this tool. The advanced argument `info-annotation-keys` is available for models trained with different sets info field annotations. In order to do this you must first train your own model with the tools CNNVariantWriteTensors and CNNVariantTrain. Otherwise, providing this argument with anything but the standard set of annotations will result in an error.

for more info visit [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037226672-CNNScoreVariants)

##### C.1 Apply a Convolutional Neural Net to filter annotated variants

```
source activate gatk
time gatk CNNScoreVariants \
--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
-I ${bam}.sorted.dup.recal.X.bam \
-V ${vcf}.X.vcf \
-R ${genome} \
-O ${vcf}.CNN.X.vcf \
--tensor-type read_tensor \
--tmp-dir ${TMP} \
   -L ${intervals} \
   -ip 100 \

```

##### C.2 Apply tranche filtering

```
gatk FilterVariantTranches \
--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
-V ${vcf}.CNN.X.vcf \
-O ${vcf}.CNN.filtered.X.vcf \
--info-key CNN_2D \
-resource ${hapmap} \
-resource ${omni} \
-resource ${G1000} \
-resource ${dbsnp} \
-resource ${mills} \
--tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \

```

#
### 3. Genotype Refinement

While every study can benefit from increased data accuracy, this workflow is especially useful for analyses that are concerned with how many copies of each variant an individual has (e.g. in the case of loss of function) or with the transmission (or de novo origin) of a variant in a family.

*pic here*
https://drive.google.com/uc?id=15fdb6LwnB4Hxw3JAUBclkJnueaV1ehV-

#### 3.1 Derive posterior probabilities of genotypes

sing the Phred-scaled genotype likelihoods (PLs) for each sample, prior probabilities for a sample taking on a HomRef, Het, or HomVar genotype are applied to derive the posterior probabilities of the sample taking on each of those genotypes. A sample’s PLs were calculated by HaplotypeCaller using only the reads for that sample. By introducing additional data like the allele counts from the 1000 Genomes project and the PLs for other individuals in the sample’s pedigree trio, those estimates of genotype likelihood can be improved based on what is known about the variation of other individuals.

SNP calls from the 1000 Genomes project capture the vast majority of variation across most human populations and can provide very strong priors in many cases. At sites where most of the 1000 Genomes samples are homozygous variant with respect to the reference genome, the probability of a sample being analyzed of also being homozygous variant is very high.

For a sample for which both parent genotypes are available, the child’s genotype can be supported or invalidated by the parents’ genotypes based on Mendel’s laws of allele transmission. Even the confidence of the parents’ genotypes can be recalibrated, such as in cases where the genotypes output by HaplotypeCaller are apparent Mendelian violations.


```
gatk \
   CalculateGenotypePosteriors \
   -V indel.SNP.recalibrated_99.9.vcf.gz \
   -ped trio_pedigree.ped \
   --supporting-callsets af-only-gnomad.hg38.vcf.gz \
   -O trio_refined_99.9.vcf.gz \
      --tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \
  
  
  #for supporting you can also use 1000genome phase 3
```

#### 3.2 Filter low quality genotypes
After the posterior probabilities are calculated for each sample at each variant site, genotypes with GQ < 20 based on the posteriors are filtered out. GQ20 is widely accepted as a good threshold for genotype accuracy, indicating that there is a 99% chance that the genotype in question is correct. Tagging those low quality genotypes indicates to researchers that these genotypes may not be suitable for downstream analysis. However, as with the VQSR, a filter tag is applied, but the data is not removed from the VCF.

```
gatk VariantFiltration \
-R ${genome} \
-V ${vcf} \
--genotype-filter-expression "GQ<20" \
--genotype-filter-name "lowGQ" \
-O ${vcf}.GQ20 \
      --tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \
```

#### 3.3 Annotate possible de novo mutations
Only If you have a trio. Dont use this for single sample

Using the posterior genotype probabilities, possible de novo mutations are tagged. Low confidence de novos have child GQ >= 10 and AC < 4 or AF < 0.1%, whichever is more stringent for the number of samples in the dataset. High confidence de novo sites have all trio sample GQs >= 20 with the same AC/AF criterion.

```
gatk VariantAnnotator \
-R reference.fasta \
-V input.vcf \
-A PossibleDeNovo \
-O output.vcf \
--tmp-dir ${TMP} \
   -L ${intervals} \
   -ip 100 \

```


## Evaluate callset <a name="Evaluate-callset"></a>


Filtering is about balancing sensitivity and precision for research aims. ​For example, genome-wide association studies can afford to maximize sensitivity over precision such that there are more false positives in the callset. Conversely, downstream analyses that require high precision, e.g. those that cannot tolerate false positive calls because validating variants is expensive, maximize precision over sensitivity such that the callset loses true positives.

visit [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035531572) for more info.

#
### Input
Variant information file (VCF) along with its index

#
### Output
Metrics table

#
### Tools
[CollectVariantCallingMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036363592-CollectVariantCallingMetrics-Picard-#--OUTPUT)  
[VariantEval](https://gatk.broadinstitute.org/hc/en-us/articles/360037224592-VariantEval-BETA-)  
[GenotypeConcordance](https://gatk.broadinstitute.org/hc/en-us/articles/360037425091-GenotypeConcordance-Picard-)  
[Hap.py](https://github.com/Illumina/hap.py)  

---
### Steps

### 1. Using hap.py
GA4GH (​Global Alliance for Genomics and Health​) recommends using ​hap.py​ for stratified variant evaluations ([1](https://docs.google.com/spreadsheets/d/1zdzNpldjYLGuuFlE_lcwDad-O7aoAUuTZkt4JmKD-Pw/edit#gid=1080374770), [2](https://docs.google.com/document/d/1jjC9TFsiDZxen0KTc2Obx6A3AHjkwAQnPV-BPhxsGn8/edit))

for more information check [Haplotype Comparison Tools](https://github.com/Illumina/hap.py) and [Hap.py User's Manual
](https://github.com/Illumina/hap.py/blob/master/doc/happy.md).

#### 1.1 Installing hap.py

Before installing hap.py, you need to install a set of requirments:

- CMake > 2.8
- GCC/G++ 4.9.2+ for compiling
- Boost 1.55+
- Python 2, version 2.7.8 or greater
- Python packages: Pandas, Numpy, Scipy, pysam, bx-python
- Java 1.8 when using vcfeval.
- boost 1.55.0
- ant 1.9.2

Download and install boost:

```
cd ~
wget http://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
tar xjf boost_1_55_0.tar.bz2
cd boost_1_55_0
./bootstrap.sh --with-libraries=filesystem,chrono,thread,iostreams,system,regex,test,program_options
./b2 --prefix=$HOME/boost_1_55_0_install install
```

install hap.py:

```
python2 install.py ~/hap.py-install --with-rtgtools

```

**GA4GH recommendation**
```
rtg vcfeval -b ${truth.vcf} -c ${query.vcf} -o /hap.py_results/rtgHC/ -t ${ref_sdf} -e ${confidence.bed} --region=${my_regions}
```
#### 1.2 Run hap.py

```
hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa
```
##### Extra options:
hese settings are the defaults in the GA4GH best-practices

```
--no-leftshift --no-decompose --engine=vcfeval
```


Restrict analysis to given (dense) regions (similar to using -T in bcftools). One example use for this is to restrict the analysis to exome-only data.

```
-T TARGETS_BEDFILE, --target-regions TARGETS_BEDFILE
```


To specify a default reference file location, you can run

```
export HGREF=path-to-your-reference.fa
```


The ```--roc``` switch specifies the feature to filter on. Hap.py translates the truth and query GQ(X) fields into the INFO fields T_GQ and Q_GQ, it tries to use GQX first, if this is not present, it will use GQ. When run without internal preprocessing any other input INFO field can be used (e.g. --roc INFO.VQSLOD for GATK).

The ```--roc-filter``` switch may be used to specify the particular VCF filter which implements a threshold on the quality score. When calculating filtered TP/FP counts, this filter will be removed, and replaced with a threshold filter on the feature specified by ```--roc```. By default, a PASS and an ALL ROC will be computed corresponding to the variant counts with all filters enabled (PASS) and no filters (ALL). When ```--roc-filter``` is specified, a third ROC curve is computed named "SEL", which shows the performance for selectively-filtered variants.

When computing precision/recall curves, we assume that higher quality scores are better, variants with scores higher than the variable threshold will "pass", all others will "fail".

The output file will be comma-separated value files giving tp / fp / fn counts, as well as precision and recall for different thresholds of the ROC feature. Here is a full example (assuming the folder hap.py contains a clone of the hap.py repository, and that hap.py can be run through PATH):
```
hap.py hap.py/example/happy/PG_NA12878_hg38-chr21.vcf.gz \
       hap.py/example/happy/NA12878-GATK3-chr21.vcf.gz \
       -f hap.py/example/happy/PG_Conf_hg38-chr21.bed.gz \
       -r hap.py/example/happy/hg38.chr21.fa \
       -o gatk-all \
       --roc QUAL --roc-filter LowQual
       
 
```

#### 1.3 Visualize
Please use the instructions given in [this page](https://github.com/Illumina/happyR).

#
### 2. Using GATK
#### 2.1. Compare callset against a known population callset
```
gatk CollectVariantCallingMetrics \
    -I filtered.vcf.gz \
    --DBSNP Homo_sapiens_assembly38.dbsnp138.vcf \
    -SD Homo_sapiens_assembly38.dict \
    -O metrics 
```
This produces detailed and summary metrics report files. The summary metrics provide cohort-level variant metrics and the detailed metrics segment variant metrics for each sample in the callset. The detail metrics give the same metrics as the summary metrics for the samples plus several additional metrics. These are explained in detail at https://broadinstitute.github.io/picard/picard-metric-definitions.html.


#### 2.2. Compare callset against a known population callset
As of this writing, VariantEval is in beta status in GATK v4.1. And so we provide an example GATK3 command, where the tool is in production status. GATK3 Dockers are available at https://hub.docker.com/r/broadinstitute/gatk3.

```
java -jar gatk3.jar \
    -T VariantEval \
    -R Homo_sapiens_assembly38.fasta \
    -eval cohort.vcf.gz \
    -D Homo_sapiens_assembly38.dbsnp138.vcf \
    -noEV \
    -EV CompOverlap -EV IndelSummary -EV TiTvVariantEvaluator \
    -EV CountVariants -EV MultiallelicSummary \
    -o cohortEval.txt
```
This produces a file containing a table for each of the evaluation modules, e.g. CompOverlap.



#### 2.3. Calculate Genotype Concordance

```
 java -jar picard.jar GenotypeConcordance \\
       CALL_VCF=input.vcf \\
       CALL_SAMPLE=sample_name \\
       O=gc_concordance.vcf \\
       TRUTH_VCF=truth_set.vcf \\
       TRUTH_SAMPLE=sample_in_truth \\
       INTERVALS=confident.interval_list \\
       MISSING_SITES_HOM_REF = true
       
```
       
       
  
## Functional annotation <a name="Functional-annotation"></a>

### Purpose

"The ultimate goal of the functional annotation process is to assign biologically relevant information to predicted polypeptides, and to the features they derive from ( e.g. gene, mRNA)." [ref](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5850084/)

#
### Input
Filtered variant information file (VCF) along with its index.

#
### Output
A table with annotations

#
### Tools
[Snpeff](pcingola.github.io/SnpEff), [Vep](https://grch37.ensembl.org/info/docs/tools/vep/index.html), [Annovar](https://annovar.openbioinformatics.org/en/latest/)

#
### Required data
A set of databases such as:

**Functional Prediction Information:**
SIFT, PolyPhen2, LRT, Mutation Taster, FATHMM, CADD & Mutation Assessor

**Disease Association:**
ClinVar, OMIM & COSMIC

**Conservation Scores:**
PhyloP, GERP++, phastCons & SiPhy

**Population Frequencies:**
1000 Genomes, GnomAD, ExAC, and Exome Variant Server

for more info [check this page](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/)

---
### Steps
Each of the mentioned softwares could be used for annotations. For more info read [this article](https://blog.goldenhelix.com/the-sate-of-variant-annotation-a-comparison-of-annovar-snpeff-and-vep/#:~:text=VEP%20and%20Annovar%20call%20this%20a%20stop%20gain%2C,signal%20for%20the%20intron%20to%20be%20spliced%20out.) although its outdated.

#
#### SnpEff
Note: this app uses databases with Ensembl chromosome name annotations. Therefore, keep in mind that you should change your chromosome names to Ensembl's version.

```
time bcftools annotate \
--rename-chrs \
ncbitoensembl.tsv \
${vcf}.CNN.filtered.dbsnp.X.vcf > ${vcf}.CNN.filtered.dbsnp.Ensembl.X.vcf
```
After this step you can annotate your vcf:

```
#snpEff
time java  -jar ~/snpEff/snpEff.jar eff \
-c /home/myassi/snpEff/snpEff.config \
-v -no-intergenic \
-i vcf \
-o vcf GRCh38.p13.RefSeq ${vcf}.CNN.filtered.dbsnp.Ensembl.X.vcf >  ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.X.vcf


#varType
time java -jar ~/snpEff/SnpSift.jar varType  ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vcf > ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.vcf
#check file
time java -jar ~/snpEff/SnpSift.jar vcfCheck  ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.X.vcf

#dbnsfp 40m
time java -jar ~/snpEff/SnpSift.jar dbnsfp -v -db ${dbnsfp} ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.X.vcf > ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.dbnsfp.X.vcf
time java -jar ~/snpEff/SnpSift.jar vcfCheck  ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.vcf

#gwas
time java -jar ~/snpEff/SnpSift.jar gwasCat -db ${gwas} ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.dbnsfp.X.vcf > ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.dbnsfp.gwas.X.vcf
```

---
#### Annovar
Before using this annotation tool, first you should prepare your data:

#
##### 1. Normalization
```
bcftools norm -m-both -o ex1.step1.vcf ex1.vcf.gz

bcftools norm -f human_g1k_v37.fasta -o ex1.step2.vcf ex1.step1.vcf
```
#
#### 2. Format Conversion

```
convert2annovar.pl -format vcf4 -includeinfo -withzyg example/ex2.vcf -outfile ex2.avinput

```
The above command takes `ex2.vcf` as input file, and generate the `ex2.avinput` as output file. The 3 extra columns are zygosity status, genotype quality and read depth.

#
#### 3. Download annotations
Before working on gene-based annotation, a gene definition file and associated FASTA file must be downloaded into a directory if they are not already downloaded. Let's call this directory as `humandb/`


##### Reference genome
```
annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene humandb/
```
Technical Notes: The above command includes -webfrom annovar, because I already pre-built the FASTA file and included them in ANNOVAR distribution site. For other gene definition systems (such as GENCODE, CCDS) or for other species (such as mouse/fly/worm/yeast), the users needs to build the FASTA file themselves. See [this page](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/) for more details.

##### Download other databases

Other databases could be download with the following code:

```
perl annotate_variation.pl --downdb --webfrom annovar --buildver hg38 [Table Name in annovar] humandb/ 
```
A list of annovar databases is provided [here](https://annovar.openbioinformatics.org/en/latest/user-guide/download/)

##### What databases to download?

**For frequency of variants in whole-exome data:** 

exac03: latest Exome Aggregation Consortium dataste with allele frequencies in ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian).

esp6500siv2: latest NHLBI-ESP project with 6500 exomes. Three separate key words are used for 3 population groupings: esp6500siv2_all, esp6500siv2_ea, esp6500siv2_aa.

gnomad_exome: allele frequency in gnomAD database whole-exome sequence data on multiple populations.

**For functional prediction of variants in whole-exome data:**

dbnsfp: this dataset already includes SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, DANN, fitCons, PhyloP and SiPhy scores, but ONLY on coding variants.

**For functional prediction of splice variants:**

dbscsnv11: dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest, which score how likely that the variant may affect splicing

spidex: deep learning based prediction of splice variants. Unlike dbscsnv11, these variants could be far away from canonical splice sites

**For disease-specific variants:**

clinvar: ClinVar database with separate columns (CLINSIG CLNDBN CLNACC CLNDSDB CLNDSDBID) for each variant (Please check the download page for the latest version, or read below for creating your own most updated version)

cosmic: the latest COSMIC database with somatic mutations from cancer and the frequency of occurence in each subtype of cancer. For more updated cosmic, see instructions below on how to make them.

icgc21: International Cancer Genome Consortium version 21 mutations.

nci60: NCI-60 human tumor cell line panel exome sequencing allele frequency data

**For variant identifiers:**

snp: dbSNP version 142

avsnp: an abbreviated version of dbSNP 142 with left-normalization by ANNOVAR developers. (Please check the download page for the latest version)

**For chromosome coordinate of each cytogenetic band:**
cytoBand

```
perl annotate_variation.pl --downdb --buildver hg38 cytoBand humandb/

```

#
#### 4. Annotate variants

Annotate variants with the `table_annovar.pl` script, which allows custom selection of multiple annotation types in the same command with specified order of the output.

Forexample:

```
perl table_annovar.pl <variant.vcf> humandb/--outfile final –buildver
hg19 --protocol refGene,cytoBand,1000g2014oct_eur,1000g2014oct_afr,exac03,
ljb26_all,clinvar_20140929,snp138 --operation g,r,f,f,f,f,f,f --vcfinput
```

`variant.vcf` refers to the name of the input VCF file.
Please follow the `--protocol` argument by the exact names of the annotation sources, and follow the `--operation` argument by the annotation type: `g` for gene-based annotation, `r` for region-based annotation and `f` for filter-based annotation. The `--outfile` argument specifies the prefix of the name of your output file.

#
#### 5. ACMG Annotation

Download and install InterVar based on the instructions [here](https://github.com/WGLab/InterVar) and [here](https://github.com/WGLab/InterVar/blob/master/docs/user-guide/manual.md)
