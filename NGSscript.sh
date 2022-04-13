#!/bin/bash

#read -p "Enter patients name: " sample
sample="10-388-MAKH"
fastq1="${sample}_1.fastq.gz"
fastq2="${sample}_2.fastq.gz"

echo making directories
mkdir -p Results/${sample}/BAM
mkdir -p Results/${sample}/calls
mkdir -p Results/${sample}/annotate

#qc
qc="Results/${sample}/qc"
#Bam
bam="Results/${sample}/BAM"
#VCF
vcf="Results/${sample}/calls"
trim="Results/${sample}/trimmed"
annotate="Results/${sample}/annotate"
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
IP=10

#PATHs
GATK="/home/myassi/RK/Downloads/gatk-4.2.6.0/gatk-package-4.2.6.0-local.jar"
export PATH="/media/myassi/01CD77B76BF0B4F0/NGS/RK/annovar.latest/annovar":$PATH
export PATH="/media/myassi/01CD77B76BF0B4F0/NGS/RK/InterVar-master":$PATH

source ~/.bashrc
source activate gatk

      echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################       Software Version         #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo
	  
echo "gatk  4.2.6.0"
conda list -n gatk | grep -E bwa\|fastqc\|multiqc\|samtools\|bcftools

      echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################            Alignment           #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

 bwa mem -M -t 4 \
   -R '@RG\tID:${sample}\tSM:${sample}\tLB:${sample}\tPU:run${sample}_1\tCN:OGENE\tPL:ILLUMINA' \
   ${genome} \
   ${paired1} \
   ${paired2} \
   | java -jar $GATK SortSam \
   -I /dev/stdin \
   -O ${bam}/${sample}.sorted.bam \
   -SO coordinate \
   --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
   --TMP_DIR ${TMP}

   [ -f ${bam}/${sample}.sorted.bam ] && echo "alignment was done succesfuly."


 samtools stats ${bam}/${sample}.sorted.bam > $qc/bamqc.txt


 multiqc Results/${sample} -f -o ${qc}


      echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################         Mark Duplicates        #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

 java -jar $GATK MarkDuplicates \
   --REMOVE_DUPLICATES false --CREATE_INDEX true \
   -I ${bam}/${sample}.sorted.bam \
   -O ${bam}/${sample}.sorted.dup.bam \
   -M ${bam}/${sample}.sorted.dup.metrics \
   --TMP_DIR ${TMP}

   [ -f ${bam}/${sample}.sorted.bam ] && echo "MarkDuplicates was done succesfuly."

   multiqc Results/${sample} -f -o ${qc}

   	  echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################        Base Recalibrator       #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

   java -jar $GATK BaseRecalibrator \
     -I ${bam}/${sample}.sorted.dup.bam \
     -R ${genome} \
     --known-sites ${dbsnp} \
     --known-sites ${mills} \
     --known-sites ${known_indels} \
     -O ${bam}/${sample}.sorted.dup.recalibration_report.table \
     -L ${intervals} \
     -ip $IP \

     [ -f ${bam}/${sample}.sorted.dup.recalibration_report.table ] && echo "BaseRecalibration was done succesfuly."

	  echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################           Apply BQSR           #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

   java -jar $GATK ApplyBQSR \
    -R ${genome} \
    -bqsr ${bam}/${sample}.sorted.dup.recalibration_report.table \
    -I ${bam}/${sample}.sorted.dup.bam  \
    -O ${bam}/${sample}.sorted.dup.recal.bam \
     -L ${intervals} \
     -ip $IP \

     multiqc Results/${sample} -f -o ${qc}

     [ -f ${bam}/${sample}.sorted.dup.recal.bam ] && echo "ApplyBQSR was done succesfuly."

	  echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################Collect Alignment Summary Metric#####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo
	 
   java -jar $GATK CollectAlignmentSummaryMetrics \
     -R ${genome} \
     -I  ${bam}/${sample}.sorted.dup.recal.bam \
     -O  ${bam}/${sample}.sorted.dup.recal.metric.alignment.tsv \
     --METRIC_ACCUMULATION_LEVEL LIBRARY

   echo percentage of reads aligned:  A good alignment if > 90%
   head -10  ${bam}/${sample}.sorted.dup.recal.metric.alignment.tsv | tail -4 | cut -f7
   #shows percentage of reads aligned
   # A good alignment if > 90%

   	  echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################        Variant Calling         #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

   time java -jar $GATK HaplotypeCaller \
   -R ${genome} \
   -I ${bam}/${sample}.sorted.dup.recal.bam \
   -O ${vcf}/${sample}.vcf \
   -bamout ${bam}/${sample}.relaligned.bam \
   --tmp-dir ${TMP} \
      -L ${intervals} \
      -ip $IP \

      [ -f ${vcf}/${sample}.vcf ] && echo "Variant Calling was done succesfuly."

	  echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################       CNN Score Variants       #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

   #CNN 84m
   time java -jar $GATK CNNScoreVariants \
   -I ${bam}/${sample}.sorted.dup.recal.bam \
   -V ${vcf}/${sample}.vcf \
   -R ${genome} \
   -O ${vcf}/${sample}.CNN.vcf \
   --tensor-type read_tensor \
   --tmp-dir ${TMP} \
      -L ${intervals} \
      -ip $IP \

      [ -f ${vcf}/${sample}.CNN.vcf ] && echo "CNNScoreVariants was done succesfuly."

      echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################    Filter Variant Tranches     #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

   #Apply tranche filtering 46m
   time java -jar $GATK FilterVariantTranches \
   -V ${vcf}/${sample}.CNN.vcf \
   -O ${vcf}/${sample}.CNN.filtered.vcf \
   --info-key CNN_2D \
   -resource ${hapmap} \
   -resource ${omni} \
   -resource ${G1000} \
   -resource ${dbsnp} \
   -resource ${mills} \
   --tmp-dir ${TMP} \
   -L ${intervals} \
      -ip $IP \

      echo
      echo "##########################################################################"
      echo "#####################                                #####################"
      echo "#####################   Filter by Genotype Quality   #####################"
      echo "#####################                                #####################"
      echo "##########################################################################"
      echo

      java -jar $GATK VariantFiltration \
      -R ${genome} \
      -V ${vcf}/${sample}.CNN.filtered.vcf \
      --genotype-filter-expression "GQ<20" \
      --genotype-filter-name "lowGQ" \
      -O ${vcf}/${sample}.CNN.filtered.GQ20.vcf \
      --tmp-dir ${TMP} \
      -L ${intervals} \
      -ip $IP \


      [ -f ${vcf}/${sample}.CNN.filtered.vcf ] && echo "CNNScoreVariants was done succesfuly."

   echo variant calling was done successfuly

   multiqc Results/${sample} -f -o ${qc}

   echo
   echo "##########################################################################"
   echo "#####################                                #####################"
   echo "#####################      Annotation Preprocess     #####################"
   echo "#####################                                #####################"
   echo "##########################################################################"
   echo

   bcftools norm -m-both \
   -o ${annotate}/${sample}.step1.vcf ${vcf}/${sample}.CNN.filtered.GQ20.vcf

   bcftools norm -f \
    ${genome} \
    -o ${annotate}/${sample}.step2.vcf \
    ${annotate}/${sample}.step1.vcf
    echo
    echo "##########################################################################"
    echo "#####################                                #####################"
    echo "#####################        Start Annotation        #####################"
    echo "#####################                                #####################"
    echo "##########################################################################"
    echo


    table_annovar.pl ${annotate}/${sample}.step2.vcf \
    annovar.latest/annovar/humandb/ -out ${annotate}/final.${sample} -buildver hg38 \
    --protocol refGene,ensGene,knownGene,rmsk,cytoBand,1000g2015aug_all,exac03,esp6500siv2_all,gnomad_genome,avsnp150,dbnsfp42a,dbscsnv11,clinvar_20210501 \
    --operation g,g,g,r,r,f,f,f,f,f,f,f,f \
    -nastring . \
    -polish \
    --remove \
    --otherinfo \
    --vcfinput

    echo
    echo "##########################################################################"
    echo "#####################                                #####################"
    echo "#####################              ACMG              #####################"
    echo "#####################                                #####################"
    echo "##########################################################################"
    echo
    Intervar.py \
      -i ${annotate}/final.${sample}.hg38_multianno.avinput \
      -t InterVar-master/intervardb \
      -o ${annotate}/final.${sample} \
      --skip_annovar -b hg38
