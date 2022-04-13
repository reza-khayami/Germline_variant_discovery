#Files--------------------
#Fastq
root="/media/myassi/01CD77B76BF0B4F0/NGS/RK/"
cd ${root}

sample="Melika_Asoodeh"
fastq1="${sample}_1.fastq.gz"
fastq2="${sample}_2.fastq.gz"

mkdir Analysis/QC/${sample}
mkdir Analysis/BAM/${sample}
mkdir Analysis/Variants/${sample}
#Trimmomatic
paired1="Analysis/QC/Melika_Asoodeh/${sample}_1.paired.fastq.gz"
unpaired1="Analysis/QC/Melika_Asoodeh/${sample}_1.upaired.fastq.gz"
paired2="Analysis/QC/Melika_Asoodeh/${sample}_2.paired.fastq.gz"
unpaired2="Analysis/QC/Melika_Asoodeh/${sample}_2.upaired.fastq.gz"

#Bam
bam="Analysis/BAM/${sample}/${sample}"



#VCF
vcf="Analysis/Variants/${sample}/${sample}"



#Datasets
genome="databases/GRCh38_index/genome_assemblies_genome_fasta/ncbi-genomes-2021-11-27/GCF_000001405.39_GRCh38.p13_genomic.fna"
dbsnp="databases/dbsnp/GCF_000001405.39.gz"
mills="databases/mills/NCBI_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_indels="databases/known_indels/NCBI_Homo_sapiens_assembly38.known_indels.vcf"
hapmap="databases/Hapmap3.3/NCBI_hapmap_3.3.hg38.vcf.gz"
G1000="databases/1000G_phase1_SNP/NCBI_1000G_phase1.snps.high_confidence.hg38.vcf.gz"
omni="databases/1000Gomni2.5/NCBI_1000G_omni2.5.hg38.vcf.gz"
dbnsfp="databases/dbnsfp/dbNSFP4.1a.txt.gz"
gwas="databases/GWAS_catalog/gwas_catalog_v1.0.2-associations_e105_r2022-01-12.tsv"
intervals="databases/Gencode/list.interval_list"
#Temp
TMP="tmp"


#Codes----------------------
#QC 4m
fastqc \
 FastQ/${fastq1} \
 FastQ/${fastq2}
 
 multiqc .
 #Trimming 8m
  java -jar ~/RK/bin/trimmomatic-0.39.jar PE \
${root}/FastQ/${fastq1} ${root}/FastQ/${fastq2} \
${paired1} \
${unpaired1} \
${paired2} \
${unpaired1} \
 ILLUMINACLIP:/home/myassi/RK/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:15 TRAILING:15 MINLEN:32 \



#QC2 4m
fastqc \
 ${paired1} \
 ${paired2}
 
#Alignment and Sorting 54m
bwa mem -M -t 2 \
  -R '@RG\tID:ZahraRahmati\tSM:ZahraRahmati\tLB:ZahraRahmati\tPU:runZahraRahmati_1\tCN:OGENE\tPL:ILLUMINA' \
  ${genome} \
  ${paired1} \
  ${paired2} \
  | java -Xmx2G -jar ~/RK/bin/picard.jar SortSam \
  -I /dev/stdin \
  -O ${bam}.sorted.bam \
  -SO coordinate \
  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 500000 \
  --TMP_DIR ${TMP}
  
  #Mark Duplicates 7m
java -Xmx2G -jar ~/RK/bin/picard.jar MarkDuplicates \
  --REMOVE_DUPLICATES false --CREATE_INDEX true \
  -I ${bam}.sorted.bam \
  -O ${bam}.sorted.dup.bam \
  -M ${bam}.sorted.dup.metrics \
  --TMP_DIR ${TMP} 
  
#make intervals

zcat gencode.v39.annotation.gtf.gz | awk '$3 == "exon" { print $1, $4, $5, $7, $18}' | tr ' ' '\t' | tr -d '";' >> interval_list2.bed
#remove the empty contigs (I did it in R)
a <- read.delim("J:/NGS/RK/databases/Gencode/interval_list2.bed", header = F)

a[rownames(a[a$V6 == "NO",] ),]
a$V6 <- ifelse(a$V2 < a$V3, "YES", "NO")

d <- a[rownames(a[a$V6 == "YES",] ),-6]

write.table(d,"J:/NGS/RK/databases/Gencode/intervals.bed", quote = F, row.names = F, col.names = F, sep = "\t")

#convert to gatk format
gatk BedToIntervalList -I interval_list.bed -O \
list.interval_list -SD GCF_000001405.39_GRCh38.p13_genomic.dict


  #Base Quality Scores Recalibration 48m + 9m +6m
 time gatk BaseRecalibrator \
   -I ${bam}.sorted.dup.bam \
   -R ${genome} \
   --known-sites ${dbsnp} \
   --known-sites ${mills} \
   --known-sites ${known_indels} \
   -O ${bam}.sorted.dup.recalibration_report.table \
   -L ${intervals} \
   -ip 100 \

gatk ApplyBQSR \
  -R ${genome} \
  -bqsr ${bam}.sorted.dup.recalibration_report.table \
  -I ${bam}.sorted.dup.bam  \
  -O ${bam}.sorted.dup.recal.bam \
   -L ${intervals} \
   -ip 100 \

gatk CollectAlignmentSummaryMetrics \
  -R ${genome} \
  -I  ${bam}.sorted.dup.recal.bam \
  -O  ${bam}.sorted.dup.recal.metric.alignment.tsv \
  --METRIC_ACCUMULATION_LEVEL LIBRARY

head -10  ${bam}.sorted.dup.recal.metric.alignment.tsv | tail -4 | cut -f7
#shows percentage of reads aligned
# A good alignment if > 90%


#Variant Calling  71m
time gatk HaplotypeCaller \
-R ${genome} \
-I ${bam}.sorted.dup.recal.bam \
-O ${vcf}.vcf.gz \
-bamout ${bam}.relaligned.bam \
--tmp-dir ${TMP} \
   -L ${intervals} \
   -ip 100 \
   

#CNN 84m
source activate gatk
time gatk CNNScoreVariants \
--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
-I ${bam}.sorted.dup.recal.bam \
-V ${vcf}.vcf.gz \
-R ${genome} \
-O ${vcf}.CNN.vcf \
--tensor-type read_tensor \
--tmp-dir ${TMP} \
   -L ${intervals} \
   -ip 100 \

#Apply tranche filtering 46m
time gatk FilterVariantTranches \
--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
-V ${vcf}.CNN.vcf \
-O ${vcf}.CNN.filtered.vcf \
--info-key CNN_2D \
-resource ${hapmap} \
-resource ${omni} \
-resource ${G1000} \
-resource ${dbsnp} \
-resource ${mills} \
--tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \

gatk CalculateGenotypePosteriors \
   -V ${vcf}.CNN.filtered.vcf \
   -O ${vcf}.CNN.filtered.post.vcf \
   -supporting /media/myassi/01CD77B76BF0B4F0/NGS/RK/databases/Gnomad/NCBI_af-only-gnomad.hg38.vcf.gz \
   --tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \
   
   
   
gatk VariantFiltration \
-R ${genome} \
-V ${vcf}.CNN.filtered.post.X.vcf \
--genotype-filter-expression "GQ<20" \
--genotype-filter-name "lowGQ" \
-O ${vcf}.CNN.filtered.post.GQ20.X.vcf \
      --tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \


 gatk VariantFiltration \
   -R ${genome} \
   -V ${vcf}.CNN.filtered.post.X.vcf \
   -O ${vcf}.CNN.filtered.post.GQ20.X.vcf \
   -G-filter "GQ >= 20" \
   --filter-name  "GQ20" \
      --tmp-dir ${TMP} \
-L ${intervals} \
   -ip 100 \
#get rsIDs 55m

time gatk  VariantAnnotator \
-R ${genome} \
--dbsnp ${dbsnp} \
-V ${vcf}.CNN.filtered.X.vcf \
-O ${vcf}.CNN.filtered.dbsnp.X.vcf \
--tmp-dir ${TMP} \
   -L ${intervals} \
   -ip 100 \
#Annotate
time bcftools annotate \
--rename-chrs \
ncbitoensembl.tsv \
${vcf}.CNN.filtered.dbsnp.X.vcf > ${vcf}.CNN.filtered.dbsnp.Ensembl.X.vcf



#snpEff 3m
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

#Evaluate 47m
#convert
time bcftools annotate \
--rename-chrs \
ensembltoncbi.tsv \
${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.vcf > ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.ncbi.vcf

time gatk CollectVariantCallingMetrics \
    -I ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.ncbi.vcf \
    --DBSNP ${dbsnp} \
    -O ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.metrics 



grep HIGH ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.vcf
java -jar ~/snpEff/SnpSift.jar extractFields -s "," -e "." ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.vcf CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].HGVS_P"


./scripts/vcfEffOnePerLine.pl \ 
java -jar ~/snpEff/SnpSift.jar extractFields -s "," -e "." ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.vcf CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT: has 'HIGH'" "ANN[*].GENE:" "ANN[*].GENEID:" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE:" "ANN[*].RANK:" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].ERRORS" "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC" > ${vcf}.new.vcf



./scripts/vcfEffOnePerLine.pl \ 
java -jar ~/snpEff/SnpSift.jar filter "ANN[*].IMPACT: has 'HIGH'" ${vcf}.CNN.filtered.dbsnp.Ensembl.snpEff.vartype.dbnsfp.gwas.vcf   > ${vcf}.new.vcf





#VQSR-------------------------------------------------------------------------------------------

#split variants
gatk SelectVariants \
-V ${vcf}.vcf \
-select-type SNP \
-O ${vcf}.vq.SNP.vcf


time gatk SelectVariants \
-V ${vcf}.vcf \
-select-type INDEL \
-O ${vcf}.vq.INDEL.vcf



time gatk VariantRecalibrator \
--java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
	-R  ${genome} \
	-V ${vcf}.vq.SNP.vcf \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	--resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap} \
	--resource:omni,known=false,training=true,truth=true,prior=12 ${omni} \
	--resource:1000G,known=false,training=true,truth=false,prior=10 ${G1000} \
	--resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp} \
	-O ${vcf}.vq.SNP.calibrate.recal \
	--tranches-file ${vcf}.vq.SNP.calibrate.cohort_snps.tranches \
	--rscript-file ${vcf}.output.plots.R



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
	--rscript-file ${vcf}.output.plots.R
