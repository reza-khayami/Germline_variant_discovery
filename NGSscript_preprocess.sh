#!/bin/bash
sample="10-388-MAKH"

#read -p "Enter patients name: " sample
fastq1="${sample}_1.fastq.gz"
fastq2="${sample}_2.fastq.gz"

echo making directories
mkdir -p Results/${sample}/qc
mkdir -p Results/${sample}/trimmed
#qc
qc="Results/${sample}/qc"
#Bam
bam="Results/${sample}/BAM"
#VCF
vcf="Results/${sample}/calls"
trim="Results/${sample}/trimmed"
#Trimmomatic
paired1="${trim}/${sample}_1.paired.fastq.gz"
unpaired1="${trim}/${sample}_1.upaired.fastq.gz"
paired2="${trim}/${sample}_2.paired.fastq.gz"
unpaired2="${trim}/${sample}_2.upaired.fastq.gz"


#Datasets
genome="databases/GRCh38_index/genome_assemblies_genome_fasta/NCBI_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
#genome="databases/GRCh38_index/genome.fa"
dbsnp="databases/dbsnp/dbsnp155ucsc.vcf.gz"
mills="databases/mills/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_indels="databases/known_indels/Homo_sapiens_assembly38.known_indels.vcf.gz"
hapmap="databases/Hapmap3.3/hapmap_3.3.hg38.vcf.gz"
G1000="databases/1000G_phase1_SNP/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
omni="databases/1000Gomni2.5/1000G_omni2.5.hg38.vcf.gz"
intervals="databases/Gencode/list.interval_list"
#Temp
TMP="tmp"

echo "###############################quality control###################################"


fastqc \
-o ${qc} \
 Samples/${fastq1} \
 Samples/${fastq2}

 multiqc ${qc} -f -o ${qc}
 google-chrome-stable ${qc}/multiqc_report.html

 read -p "Do you want to start Trimmomatic? (y/n)" -n 1 -r
 echo    # (optional) move to a new line
 if [[ ! $REPLY =~ ^[Yy]$ ]]
 then
     exit 1
 fi
echo "###############################start trimming###################################"

 java -jar ~/RK/bin/trimmomatic-0.39.jar PE \
Samples/${fastq1} Samples/${fastq2} \
 ${paired1} \
 ${unpaired1} \
 ${paired2} \
 ${unpaired1} \
 ILLUMINACLIP:/home/myassi/RK/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:15 TRAILING:15 MINLEN:32

 echo "###############################post trim quality control###################################"

 #QC2 4m
 fastqc \
 -o ${qc} \
 ${paired1} \
 ${paired2}

 multiqc ${qc} -f -o ${qc}

 [ -f "${paired1}" -a -f "${paired2}" ] && echo "Preprocessing was done succesfuly."


google-chrome-stable ${qc}/multiqc_report.html
