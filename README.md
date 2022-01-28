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
 
![](https://doc-00-90-docs.googleusercontent.com/docs/securesc/ukbe88t0fjafimk2ma5n0e2lm57nhhpt/7ngug6d52paulu5flu5hqoq6gs420cu6/1643391600000/03300952448340018413/10869494293900274643/1HKtzOeobgOVjCXEUE0-5378ocBz6Age7?authuser=0&nonce=53ichpq6b40se&user=10869494293900274643&hash=1jukoceo9ivoqevjb9hs72boalnd5dfn)



## Data preprocessing <a name="Data-preprocessing"></a>
### Purpose
The is the obligatory first phase that must precede all variant discovery. It involves pre-processing the raw sequence data (provided in FASTQ format) to produce analysis-ready BAM files. This involves alignment to a reference genome as well as some data cleanup operations to correct for technical biases and make the data suitable for analysis.

### Input
FASTQ format, Compressed or uncompressed, single end or paired end.

### Output
BAM file and its index file ready for variant discovery

### Tools
BWA, SAMtools, Bcftools, GATK, multiqc


### References
-----------------------------------------------------
### Steps
#### 1. Quality control <a name="Quality-control"></a>
#### 2. Trimming <a name="Trimming"></a>
#### 3. Alignment <a name="Alignment"></a>
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

