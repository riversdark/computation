#!/bin/bash

# Script to call germline and somatic variants in a human WGS paired end reads 2 X 100bp
# germline variant discovery - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# somatic variant discovery - https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-

######################################################################
# 0. prepare reference files
echo "Step 0: prepare reference files..."

REF_DIR="/mnt/nas/wgs/hg38"
REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# download reference files if not already present
if [ ! -f ${REF_DIR}/hg38.fa ]; then
    echo "Downloading reference files..."
    wget -P ${REF_DIR} ${REF_URL}
    gunzip ${REF_DIR}/hg38.fa.gz
fi

# index reference files if not already present
if [ ! -f ${REF_DIR}/hg38.fa.fai ]; then
    echo "Indexing reference files..."
    samtools faidx ${REF_DIR}/hg38.fa
fi

# BWA index reference files if not already present
if [ ! -f ${REF_DIR}/hg38.fa.amb ]; then
    echo "BWA indexing reference files..."
    bwa index ${REF_DIR}/hg38.fa
fi

# create sequence dictionary if not already present
if [ ! -f ${REF_DIR}/hg38.dict ]; then
    echo "Creating sequence dictionary..."
    gatk CreateSequenceDictionary -R ${REF_DIR}/hg38.fa -O ${REF_DIR}/hg38.dict
fi

# download known sites files for BQSR from GATK resource bundle if not already present
KNOWN_SITES_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_SITES_IDX_URL="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
KNOWN_SITES="Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_SITES_IDX="Homo_sapiens_assembly38.dbsnp138.vcf.idx"
if [ ! -f ${REF_DIR}/${KNOWN_SITES} ]; then
    echo "Downloading known sites files..."
    wget -P ${REF_DIR} ${KNOWN_SITES_URL}
    wget -P ${REF_DIR} ${KNOWN_SITES_IDX_URL}
fi

# list the reference files
echo "Reference files:"
ls -hlt ${REF_DIR}
echo "Reference files preparation finished."

# ######################################################################
# 1. Pre-processing of raw sequencing data
echo "Step 1: Pre-processing of raw sequencing data..."

READS_DIR="/mnt/nas/wgs/reads/SRR062634"
READS_FILE1="SRR062634_1.filt"
READS_FILE2="SRR062634_2.filt"

# a. Quality control of raw sequencing data using FastQC, if not already present
if [ ! -f ${READS_DIR}/${READS_FILE1}_fastqc.zip ]; then
    echo "Running FastQC on ${READS_FILE1}..."
    fastqc ${READS_DIR}/${READS_FILE1}.fastq.gz -o ${READS_DIR}
fi

if [ ! -f ${READS_DIR}/${READS_FILE2}_fastqc.zip ]; then
    echo "Running FastQC on ${READS_FILE2}..."
    fastqc ${READS_DIR}/${READS_FILE2}.fastq.gz -o ${READS_DIR}
fi

echo "FastQC finished."

# b. Adapter trimming using Trimmomatic


# c. Alignment to a reference genome using BWA-MEM, if not already present
if [ ! -f ${READS_DIR}/${READS_FILE1}.sam ]; then
    echo "Running BWA-MEM on ${READS_FILE1}.fastq.gz..."
    bwa mem -t 50 -M -R "@RG\tID:SRR062634\tSM:SRR062634\tPL:ILLUMINA" ${REF_DIR}/hg38.fa ${READS_DIR}/${READS_FILE1}.fastq.gz ${READS_DIR}/${READS_FILE2}.fastq.gz > ${READS_DIR}/${READS_FILE1}.sam
fi

# Show the flag statistics
echo "Flag statistics:"
# samtools flagstat ${READS_DIR}/${READS_FILE1}.sam

# ######################################################################
# 2. Sort and index the aligned BAM files.
echo "Step 2: Sort and index the aligned BAM files..."

# a. Sort the aligned BAM files
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted.bam ]; then
    echo "Sorting ${READS_FILE1}.sam..."
    samtools sort -@ 50 -o ${READS_DIR}/${READS_FILE1}_sorted.bam ${READS_DIR}/${READS_FILE1}.sam
fi

# b. Index the sorted BAM files
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted.bam.bai ]; then
    echo "Indexing ${READS_FILE1}_sorted.bam..."
    samtools index -@ 50 ${READS_DIR}/${READS_FILE1}_sorted.bam
fi

echo "Sorting and indexing finished."

# ######################################################################
# 3. Mark duplicates to identify and flag any potential PCR duplicates.
echo "STEP 3: Mark Duplicates"

# a. Mark duplicates
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup.bam ]; then
    echo "Marking duplicates in ${READS_FILE1}_sorted.bam..."
    gatk MarkDuplicates -I ${READS_DIR}/${READS_FILE1}_sorted.bam -O ${READS_DIR}/${READS_FILE1}_sorted_dedup.bam --METRICS_FILE ${READS_DIR}/${READS_FILE1}_sorted_dedup.metrics.txt
fi

# b. Index the marked duplicates BAM files
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup.bam.bai ]; then
    echo "Indexing ${READS_FILE1}_sorted_dedup.bam..."
    samtools index -@ 50 ${READS_DIR}/${READS_FILE1}_sorted_dedup.bam
fi

echo "Marking duplicates finished."


# ######################################################################
# 4. Base quality score recalibration to correct for systematic biases in the sequencing data.
echo "STEP 4: Base quality recalibration"

# a. build the model
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_data.table ]; then
    echo "Building the model..."
    gatk BaseRecalibrator -R ${REF_DIR}/hg38.fa -I ${READS_DIR}/${READS_FILE1}_sorted_dedup.bam -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_data.table --known-sites ${REF_DIR}/${KNOWN_SITES}
fi


# b. Apply the model to adjust the base quality scores
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal.bam ]; then
    echo "Applying the model..."
    gatk ApplyBQSR -R ${REF_DIR}/hg38.fa -I ${READS_DIR}/${READS_FILE1}_sorted_dedup.bam -bqsr ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_data.table -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal.bam
fi

echo "Base quality recalibration finished."


# ######################################################################
# # 5. Collect alignment and insert size metrics to assess the quality of the alignment and to detect potential problems with the library preparation.
echo "STEP 5: Collect Alignment & Insert Size Metrics"

# a. Collect alignment metrics
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_alignment_metrics.txt ]; then
    echo "Collecting alignment metrics..."
    gatk CollectAlignmentSummaryMetrics -R ${REF_DIR}/hg38.fa -I ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal.bam -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_alignment_metrics.txt -H ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_alignment_histogram.pdf
fi

# b. Collect insert size metrics
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_insert_size_metrics.txt ]; then
    echo "Collecting insert size metrics..."
    gatk CollectInsertSizeMetrics -R ${REF_DIR}/hg38.fa -I ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal.bam -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_insert_size_metrics.txt -H ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_insert_size_histogram.pdf
fi

echo "Collecting alignment & insert size metrics finished."


# ######################################################################
# 6. Calling variants with HaplotypeCaller to generate a raw variant call set in VCF format.
echo "STEP 6: Call Variants - gatk haplotype caller"

# a. Call variants
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_variants.vcf ]; then
    echo "Calling variants..."
    gatk HaplotypeCaller -R ${REF_DIR}/hg38.fa -I ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal.bam -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_variants.vcf
fi

echo "Calling variants finished."

# ######################################################################
# 7. Extract SNPs and INDELs from the raw variant call set.
echo "STEP 7: Extract SNPs and INDELs"

# a. Extract SNPs
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_snps.vcf ]; then
    echo "Extracting SNPs..."
    gatk SelectVariants -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_variants.vcf -select-type SNP -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_snps.vcf
fi

# b. Extract INDELs
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_indels.vcf ]; then
    echo "Extracting INDELs..."
    gatk SelectVariants -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_variants.vcf -select-type INDEL -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_indels.vcf
fi

echo "Extracting SNPs and INDELs finished."


# ######################################################################
# 8. Filter SNPs and INDELs to remove low quality calls.
echo "STEP 8: Filter SNPs and INDELs"

# a. Filter SNPs
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps.vcf ]; then
    echo "Filtering SNPs..."
    gatk VariantFiltration -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_snps.vcf -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "my_snp_filter"
fi

# b. Filter INDELs
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_indels.vcf ]; then
    echo "Filtering INDELs..."
    gatk VariantFiltration -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_raw_indels.vcf -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_indels.vcf --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "my_indel_filter"
fi

echo "Filtering SNPs and INDELs finished."

# ######################################################################
# 9. Merge SNPs and INDELs into a single VCF file.
echo "STEP 9: Merge SNPs and INDELs"

# a. Merge SNPs and INDELs
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels.vcf ]; then
    echo "Merging SNPs and INDELs..."
    gatk MergeVcfs -I ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps.vcf -I ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_indels.vcf -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels.vcf
fi

echo "Merging SNPs and INDELs finished."

# ######################################################################
# 10. Annotate variants with Funcotator
echo "STEP 10: Annotate variants with Funcotator"

DATA_SOURCES_VERSION="v1.7.20200521g"
# a. Download data sources for variant annotation
if [ ! -d ${REF_DIR}/funcotator_dataSources.${DATA_SOURCES_VERSION} ]; then
    echo "Downloading data sources for variant annotation..."
    gatk FuncotatorDataSourceDownloader  --germline --validate-integrity --extract-after-download -O ${REF_DIR}/funcotator_dataSources.${DATA_SOURCES_VERSION}
    # gatk FuncotatorDataSourceDownloader  --somatic --validate-integrity --extract-after-download -O ${REF_DIR}/funcotator_dataSources
fi

# b. Annotate variants
if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated.vcf ]; then
    echo "Annotating variants..."
    gatk Funcotator -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels.vcf -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated.vcf --data-sources-path ${REF_DIR}/funcotator_dataSources.${DATA_SOURCES_VERSION} --ref-version hg38 --output-file-format VCF
fi

echo "Annotating variants finished."

# ######################################################################
# 11. Create a VCF file with only the PASS variants
echo "STEP 11: Create a VCF file with only the PASS variants"

if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated_PASS.vcf ]; then
    echo "Creating a VCF file with only the PASS variants..."
    gatk SelectVariants -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated.vcf -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated_PASS.vcf --exclude-filtered
fi

echo "Creating a VCF file with only the PASS variants finished."

# ######################################################################
# 12. Convert VCF file to TSV file
echo "STEP 12: Convert VCF file to TSV file"

if [ ! -f ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated_PASS.tsv ]; then
    echo "Converting VCF file to TSV file..."
    gatk VariantsToTable -R ${REF_DIR}/hg38.fa -V ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated_PASS.vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F DP -F AF -F AC -F AN -F BaseQRankSum -F ClippingRankSum -F DP -F ExcessHet -F FS -F InbreedingCoeff -F MQ -F MQRankSum -F QD -F ReadPosRankSum -F SOR -F Funcotator -O ${READS_DIR}/${READS_FILE1}_sorted_dedup_recal_filtered_snps_indels_funcotated_PASS.tsv
fi

echo "Converting VCF file to TSV file finished."

# ######################################################################
# ######################################################################
# ######################################################################
echo "All done!"

ls -hlt ${READS_DIR}

