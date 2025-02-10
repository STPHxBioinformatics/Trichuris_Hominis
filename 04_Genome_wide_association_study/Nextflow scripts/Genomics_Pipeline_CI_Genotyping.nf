nextflow.enable.dsl=2
params.reads = "/scicore/home/schpie00/GROUP/20230601_Raw_sequencing_Trichuris_whole_genome/2023*/GFB-*TRI*_R{1,2}_001.fastq.gz"
params.trimmomatic = "/scicore/home/schpie00/baer0006/IluminiPilot/Nextflow/packages/Trimmomatic-0.39/trimmomatic-0.39.jar"
params.outdir = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/Nextflow/"
params.reference = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/Reference/Trichuris_cote_divoire_freeze.fasta"
params.reference_path = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/Reference/Trichuris_cote_divoire_freeze"
params.path_kraken2_db = "/scicore/data/managed/.store/kraken_191029/"
params.flagstat_output = "/scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/flagstat/"


log.info """\
	
	=============================================
	SWISS-TPH WHOLE-GENOME-PIPELINE: CAN'T RESIST
	=============================================
	reads directory: 	        ${params.reads}
	reference file location:    ${params.reference}
	"""
	.stripIndent()


process trimming {
	
	input:
	tuple val(pair_id), path(reads)
	
	output:
	tuple val(pair_id), path("${pair_id}_trimmed_R{1,2}.fastq")
	
	script:
	"""
	java -jar $baseDir/packages/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33\
		${reads[0]}\
		${reads[1]}\
		${pair_id}_trimmed_R1.fastq\
		${pair_id}_trimmed_R1_unpaired.fastq\
		${pair_id}_trimmed_R2.fastq\
		${pair_id}_trimmed_R2_unpaired.fastq\
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
	"""
}

process index {
	
	input:
	path(reference)

	output:
	val(reference.baseName)

	script:
	"""
	ml load BWA/0.7.17-GCC-10.3.0
	ml load SAMtools/1.14-GCC-10.3.0
	bwa index ${reference};
	
	samtools faidx ${reference};
	samtools dict ${reference} > ${params.reference_path}.dict
	"""
}


process alignment {
	maxForks 60

	input:
	val(reference)
	tuple val(pair_id), path(reads)
	
	output:
	tuple val(pair_id), path("${pair_id}.bam")
	
	script:
	readGroup = "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:ILLUMINA\\tPM:NOVASEQ\\tSM:${pair_id}"
	"""
	ml load BWA/0.7.17-GCC-10.3.0
	ml load SAMtools/1.14-GCC-10.3.0
	bwa mem \
	-t 4 $params.reference \
	-R \"${readGroup}\" \
	${reads[0]} \
	${reads[1]} | samtools view --threads 4 -b - \
	> ${pair_id}.bam
	"""
}


process sorting {
	maxForks 20

	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}_filtered.bam")

	script:
	"""
	ml load SAMtools/1.14-GCC-10.3.0
	samtools sort -o ${pair_id}_sorted.bam $reads;
	samtools view -F 4 -b -o ${pair_id}_filtered.bam ${pair_id}_sorted.bam;
	samtools index ${pair_id}_filtered.bam ${pair_id}_filtered.bam.bai
	"""
}

process flagstat {
	maxForks 20
	
	input:
	tuple val(pair_id), path(reads)

	output: 
	file("${pair_id}.flagstat")

	script:
	"""
	ml load SAMtools/1.14-GCC-10.3.0
	samtools flagstat ${reads} > ${pair_id}.flagstat
	"""
}

process multiqc_flagstat {
	publishDir "$params.flagstat_outputrun2", mode: 'copy'

	input:
	file(flagstat)

	output:
	file "multiqc_report_flagstat.html"
	
	script:
	"""
	ml load MultiQC/1.11-foss-2018b-Python-3.6.6
	multiqc .;
	mv multiqc_report.html multiqc_report_flagstat.html
	"""
}

process kraken2_impurities {
	maxForks 2
	
	input:
	tuple val(pair_id), path(reads)
	
	output:
	file("${pair_id}_kraken2report")
	
	script:
	"""
	ml load Kraken2/2.1.1-foss-2018b-Perl-5.28.0
	kraken2 \
	--db $params.path_kraken2_db \
	--report ${pair_id}_kraken2report \
	--paired ${reads[0]} ${reads[1]}
	"""
}

process multiqc_kraken2 {
	publishDir "$params.flagstat_output*", mode: 'copy'

	input:
	file(kraken2)

	output:
	file "multiqc_report_kraken2.html"
	
	script:
	"""
	ml load MultiQC/1.11-foss-2018b-Python-3.6.6
	multiqc .;
	mv multiqc_report.html multiqc_report_kraken2.html
	"""
}

process mark_duplicates_spark {
	maxForks 5
	
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}_sorted_dedup.bam")
	
	script:
	"""
	ml purge
	ml load Java/8u212b03
	ml load GATK/4.2.6.1-foss-2018b-Java-1.8
	gatk MarkDuplicatesSpark \
	-I ${reads} \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam 
	"""
}

process classify_gender {
    input:
    tuple val(pair_id), path(bam_file)

    output:
    path("${pair_id}.txt")

    script:
    """
    ml load SAMtools/1.16.1-GCC-10.3.0
    samtools index ${bam_file};
    cov_1="\$(samtools coverage -r TTRE_chr1_scaffold1 ${bam_file} | awk 'NR==2 {print \$7}')";
    cov_2="\$(samtools coverage -r TTRE_chr2_scaffold39 ${bam_file} | awk 'NR==2 {print \$7}')";
    ratio=\$(echo "\$cov_1 / \$cov_2" | bc -l)
    echo "\$cov_1\n\$cov_2\n\$ratio" >> ${pair_id}.txt
    """

}

process merge_gender_data {
    input:
    path(gender_files)

    output:
    path("merged_coverage.txt")

    script:
    """
    ml load Python/3.9.5-GCCcore-10.3.0
    python /scicore/home/schpie00/baer0006/Trichuris_Trichiura_CI/analysis_pipeline/merge_coverage.py
    """

}

process initial_basecalling_BQSR {
	maxForks 56
	
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}_first_raw_variants.vcf")
	
	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11 
	gatk HaplotypeCaller \
	-R ${params.reference} \
	-I ${reads} \
	--heterozygosity 0.015 \
    	--indel-heterozygosity 0.01 \
	-O ${pair_id}_first_raw_variants.vcf
	"""
}

process initial_SNP_BQSR {
	maxForks 56
	
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}_SNP_FOR_BQSR.vcf"), path("${pair_id}_SNP_FOR_BQSR.vcf.idx")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	
	gatk SelectVariants \
	-R ${params.reference} \
	-V ${reads} \
	--select-type-to-include SNP \
	-O ${pair_id}_first_raw_variants_SNP.vcf;

	gatk VariantFiltration \
   	-R ${params.reference} \
   	-V ${pair_id}_first_raw_variants_SNP.vcf \
	--filter-expression ' QUAL < 194 || QD < 10.3 || MQ < 43.61 || FS > 9.9 || SOR > 2.30 || MQRankSum < -1.269 || ReadPosRankSum < -1.174 ' \
	--filter-name "SNP_5_filtered" \
	-O ${pair_id}_filtered_first_raw_variants_SNP.vcf;

	gatk SelectVariants \
	--exclude-filtered \
	-V ${pair_id}_filtered_first_raw_variants_SNP.vcf \
	-O ${pair_id}_SNP_FOR_BQSR.vcf
	"""
}

process initial_INDEL_BQSR {
	maxForks 56
	
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("${pair_id}_INDEL_FOR_BQSR.vcf"), path("${pair_id}_INDEL_FOR_BQSR.vcf.idx")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	gatk SelectVariants \
	-R ${params.reference} \
	-V ${reads} \
	--select-type-to-include INDEL \
	-O ${pair_id}_first_raw_variants_INDEL.vcf;

	gatk VariantFiltration \
   	-R ${params.reference} \
   	-V ${pair_id}_first_raw_variants_INDEL.vcf \
	--filter-expression ' QUAL < 101 || QD < 5.3 || FS > 11.4 || ReadPosRankSum < -1.336 ' \
	--filter-name "INDEL_5_filtered" \
	-O ${pair_id}_filtered_first_raw_variants_INDEL.vcf;
	
	gatk SelectVariants \
	--exclude-filtered \
	-V ${pair_id}_filtered_first_raw_variants_INDEL.vcf \
	-O ${pair_id}_INDEL_FOR_BQSR.vcf
	"""
}

process BQSR {
	maxForks 56

	input:
	tuple val(pair_id), path(SNPs), path(SNPindex), path(INDELs), path(INDELindex), path(BAM)

	output:
	tuple val(pair_id), path("${pair_id}_BQSR_recal.bam"), path("${pair_id}_BQSR_data_recalibration_plots.pdf"), path("${pair_id}_BQSR_data.table"), path("${pair_id}_BQSR_data_post.table")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	ml load  R/4.2.2-foss-2021a
	
	gatk BaseRecalibrator \
	-R ${params.reference} \
	-I ${BAM} \
	--known-sites ${SNPs} \
	--known-sites ${INDELs} \
	-O ${pair_id}_BQSR_data.table;

	gatk ApplyBQSR \
	-R ${params.reference} \
	-I ${BAM} \
	-bqsr ${pair_id}_BQSR_data.table\
	-O ${pair_id}_BQSR_recal.bam;

	gatk BaseRecalibrator \
	-R ${params.reference} \
	-I ${pair_id}_BQSR_recal.bam \
	--known-sites ${SNPs} \
	--known-sites ${INDELs} \
	-O ${pair_id}_BQSR_data_post.table;

	gatk AnalyzeCovariates \
	-before ${pair_id}_BQSR_data.table \
	-after ${pair_id}_BQSR_data_post.table \
	-plots ${pair_id}_BQSR_data_recalibration_plots.pdf
	"""
}

process GVCF {
	maxForks 56

	input:
	tuple val(pair_id), path(recal_BAMs), path(recalplot), path(datatable_pre), path(datatable_post)

	output:
	tuple path("${pair_id}.g.vcf.gz"), path("${pair_id}.g.vcf.gz.tbi")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	
	gatk HaplotypeCaller \
	-R ${params.reference} \
	-I ${recal_BAMs} \
	--heterozygosity 0.015 \
   	--indel-heterozygosity 0.01 \
	-O ${pair_id}.g.vcf.gz \
	-A DepthPerAlleleBySample \
	-A Coverage \
	-A ExcessHet \
	-A FisherStrand \
	-A MappingQualityRankSumTest \
	-A StrandOddsRatio \
	-A RMSMappingQuality \
	-A ReadPosRankSumTest \
	-A DepthPerSampleHC \
	-A QualByDepth \
	-ERC GVCF
	"""
}

process merge_GVCF {
	input:
	path(VCFs)
	path(VCF_indices)

	output:
	tuple path("cohort.g.vcf.gz"), path("cohort.g.vcf.gz.tbi")

	script:
	
	def input_list = VCFs.collect{"--variant $it"}.join(' ')
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11

	gatk CombineGVCFs \
	-R ${params.reference} \
	$input_list \
	-O cohort.g.vcf.gz
	"""
}

process genotyping {
	input:
	tuple path(gVCF_cohort), path(index)

	output:
	tuple path("Trichuris_CI_cohort.vcf.gz"), path("Trichuris_CI_cohort.vcf.gz.tbi")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11

	gatk GenotypeGVCFs \
	-R ${params.reference} \
	-V ${gVCF_cohort} \
	-A DepthPerAlleleBySample \
	-A Coverage \
	-A ExcessHet \
	-A FisherStrand \
	-A MappingQualityRankSumTest \
	-A StrandOddsRatio \
	-A RMSMappingQuality \
	-A ReadPosRankSumTest \
	-A DepthPerSampleHC \
	-A QualByDepth \
	--heterozygosity 0.015 \
    	--indel-heterozygosity 0.01 \
	-O Trichuris_CI_cohort.vcf.gz
	"""
}

process select_geno_SNP_INDEL {
	input:
	tuple path(VCF_cohort), path(index)
	
	output:
	tuple path("Trichuris_CI_cohort.genoSNP.vcf.gz"), path("Trichuris_CI_cohort.genoSNP.vcf.gz.tbi"), path("Trichuris_CI_cohort.genoINDEL.vcf.gz"), path("Trichuris_CI_cohort.genoINDEL.vcf.gz.tbi")
	
	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	gatk SelectVariants \
	-R ${params.reference} \
	--variant ${VCF_cohort} \
	--select-type-to-include SNP \
	-O Trichuris_CI_cohort.genoSNP.vcf.gz
	
	gatk SelectVariants \
	-R ${params.reference} \
	--variant ${VCF_cohort} \
	--select-type-to-include INDEL \
	-O Trichuris_CI_cohort.genoINDEL.vcf.gz
	"""
}

process make_table_geno_SNP_INDEL {
	input:
	tuple path(genoSNP_vcf), path(genoSNP_vcf_index), path(genoINDEL_vcf), path(genoINDEL_vcf_index)

	output:
	tuple path("geno_snp_vcf.table"), path("geno_indel_vcf.table")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	gatk VariantsToTable \
	-R ${params.reference} \
	--variant ${genoSNP_vcf} \
	--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
	-O geno_snp_vcf.table


	gatk VariantsToTable \
	-R ${params.reference} \
	--variant ${genoINDEL_vcf} \
	--fields CHROM --fields POS --fields QUAL --fields QD --fields DP --fields MQ --fields MQRankSum --fields FS --fields ReadPosRankSum --fields SOR \
	-O geno_indel_vcf.table
	"""
}

process filter_geno_SNP_INDEL {
	input:
	tuple path(genoSNP_vcf), path(genoSNP_vcf_index), path(genoINDEL_vcf), path(genoINDEL_vcf_index)

	output:
	path ("Trichuris_CI_cohort.geno.final.recode.vcf")

	script:
	"""
	ml load  GATK/4.2.4.1-GCCcore-10.3.0-Java-11
	ml load VCFtools/0.1.16-GCC-10.3.0

	gatk VariantFiltration \
	-R ${params.reference} \
	--variant ${genoSNP_vcf} \
	--filter-expression ' QUAL < 35 || QD < 1.3 || MQ < 28.76 || FS > 50.5 || SOR > 2.30 || MQRankSum < -4.261 || ReadPosRankSum < -1.421 ' \
	--filter-name "SNP_filtered" \
	-O Trichuris_CI_cohort.genoSNP.filtered.vcf
	
	gatk VariantFiltration \
	-R ${params.reference} \
	--variant ${genoINDEL_vcf} \
	--filter-expression ' QUAL < 32 || QD < 0.5 || FS > 30.1 || ReadPosRankSum < -1.836 ' \
	--filter-name "INDEL_filtered" \
	-O Trichuris_CI_cohort.genoINDEL.filtered.vcf;
	
	gatk MergeVcfs \
	-I Trichuris_CI_cohort.genoSNP.filtered.vcf \
	-I Trichuris_CI_cohort.genoINDEL.filtered.vcf \
	-O Trichuris_CI_cohort.genoMERGED.filtered.vcf;
	
	gatk VariantFiltration \
	-R ${params.reference} \
	--variant Trichuris_CI_cohort.genoMERGED.filtered.vcf \
	--genotype-filter-expression ' DP < 3 '  \
	--genotype-filter-name "DP_genotype_3" \
	--output Trichuris_CI_cohort.geno.genotype_filtered.vcf;
	
	gatk SelectVariants \
	-R ${params.reference} \
	--variant Trichuris_CI_cohort.geno.genotype_filtered.vcf \
	--set-filtered-gt-to-nocall \
	--output Trichuris_CI_cohort.geno.genotype_filtered_nocall.vcf;
	
	vcftools \
	--vcf Trichuris_CI_cohort.geno.genotype_filtered_nocall.vcf \
	--remove-filtered-geno-all \
	--remove-filtered-all \
	--min-alleles 2 \
	--max-alleles 2 \
	--maf 0.02 \
	--recode \
	--recode-INFO-all \
	--out Trichuris_CI_cohort.geno.final
	
	vcftools \
	--vcf Trichuris_CI_cohort.geno.final.recode.vcf \
	--remove-indels
	"""
}


workflow {
	read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)

	index(params.reference)
	index_ch = index.out

	alignment(index_ch, read_pairs_ch)
	aligned_ch = alignment.out

	flagstat(aligned_ch)
	flagstat_ch = flagstat.out

	multiqc_flagstat(flagstat_ch.collect())

	sorting(aligned_ch)
	sorted_filtered_ch = sorting.out
	
	mark_duplicates_spark(sorted_filtered_ch)
	deduplicated_reads_ch = mark_duplicates_spark.out

	classify_gender(sorted_filtered_ch)
	classify_gender_ch = classify_gender.out
	
	merge_gender_data(classify_gender_ch.collect())
	merge_gender_data_ch = merge_gender_data.out
	
	initial_basecalling_BQSR(deduplicated_reads_ch)
	first_raw_variants_ch = initial_basecalling_BQSR.out
	
	initial_SNP_BQSR(first_raw_variants_ch)
	SNPs_for_BQSR_ch = initial_SNP_BQSR.out

	initial_INDEL_BQSR(first_raw_variants_ch)
	INDELs_for_BQSR_ch = initial_INDEL_BQSR.out

	INDEL_SNP_BQSR_ch = SNPs_for_BQSR_ch.join(INDELs_for_BQSR_ch).join(deduplicated_reads_ch)

	BQSR(INDEL_SNP_BQSR_ch)
	BQSR_bam_ch = BQSR.out

	GVCF(BQSR_bam_ch)
	GVCF_ch = GVCF.out
	
	merge_GVCF(GVCF_ch.map { it[0] }.collect(), GVCF_ch.map { it[1] }.collect())
	merged_GVCF_ch = merge_GVCF.out

	genotyping(merged_GVCF_ch)
	genotyped_merged_vcf_ch = genotyping.out

    select_geno_SNP_INDEL(genotyped_merged_vcf_ch)
	geno_SNP_INDEL_ch = select_geno_SNP_INDEL.out
	
	 make_table_geno_SNP_INDEL(geno_SNP_INDEL_ch)
	 make_table_geno_SNP_INDEL_ch = make_table_geno_SNP_INDEL.out

	 filter_geno_SNP_INDEL(geno_SNP_INDEL_ch)
	 filter_geno_SNP_INDEL_ch = filter_geno_SNP_INDEL.out
}