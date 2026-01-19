#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Input parameters (must be specified via config file or command line)
params.project = "cohort"
params.bam_glob = null
params.out_dir = "out_dir"
params.ref = null
params.snp_site_vcf = null
params.indel_site_vcf = null
params.snp_markers_intervals = null
params.indel_markers_intervals = null
params.sv_sites_vcf = null

// Optional tool paths - if not specified, tools will be used from system PATH
params.gatk_java_path = null  // Java 8 path for GATK UnifiedGenotyper only (required for GATK 3.7)
params.java_path = null       // Java path for other tools (Picard, Beagle, etc.)
params.delly_path = null
params.bcftools_path = null
params.samtools_path = null
params.tabix_path = null

// Required software paths (JAR files)
params.gatk_path = "./GenomeAnalysisTK3.7.jar"
params.picard_path = null
params.beagle_path = null

// Resource configuration (with defaults)
params.threads = 48
params.gatk_memory = '100g'
params.beagle_memory = '300 GB'
params.beagle_cpus = 64

// Fixed path for assign_id script (always in ./scripts/assign_id.sh)
def assign_id_path = './scripts/assign_id.sh'
params.assign_id = file(assign_id_path).toAbsolutePath()

// Convert relative paths to absolute paths (if specified)
if (params.gatk_path) {
    params.gatk = file(params.gatk_path).toAbsolutePath()
} else {
    params.gatk = null
}
if (params.beagle_path) {
    params.beagle = file(params.beagle_path).toAbsolutePath()
} else {
    params.beagle = null
}
if (params.picard_path) {
    params.picard = file(params.picard_path).toAbsolutePath()
} else {
    params.picard = null
}

// Set tool executables - use specified path if provided, otherwise use tool name (from PATH)
// Handle both null and empty string cases
// GATK Java (Java 8) - for GATK UnifiedGenotyper only (GATK 3.7 requires Java 8)
params.gatk_java = (params.gatk_java_path && params.gatk_java_path.toString().trim()) ? file(params.gatk_java_path).toAbsolutePath().toString() : 'java'
// General Java - for Picard, Beagle, and other tools (can use newer Java versions)
params.java = (params.java_path && params.java_path.toString().trim()) ? file(params.java_path).toAbsolutePath().toString() : 'java'
params.delly = (params.delly_path && params.delly_path.toString().trim()) ? file(params.delly_path).toAbsolutePath().toString() : 'delly'
params.bcftools = (params.bcftools_path && params.bcftools_path.toString().trim()) ? file(params.bcftools_path).toAbsolutePath().toString() : 'bcftools'
params.samtools = (params.samtools_path && params.samtools_path.toString().trim()) ? file(params.samtools_path).toAbsolutePath().toString() : 'samtools'
params.tabix = (params.tabix_path && params.tabix_path.toString().trim()) ? file(params.tabix_path).toAbsolutePath().toString() : 'tabix'

include { INDEX_REFERENCE; UNIFIED_GENOTYPER_SNP; UNIFIED_GENOTYPER_INDEL; GATK_SNP_FORMAT; GATK_INDEL_FORMAT } from './modules/snp_indel_gt'
include { DELLY_SV_GENOTYPE; BCFTOOLS_MERGE_GENOTYPE } from './modules/sv_gt'
include { CONCAT_VCF; BEAGLE_IMPUTATION; POP_SNP; POP_INDEL; POP_SV } from './modules/utils'

// Validate required parameters
if( !params.project ) error "params.project must be specified"
if( !params.bam_glob ) error "params.bam_glob (BAM file glob pattern) must be specified"
if( !params.out_dir ) error "params.out_dir (output directory) must be specified"
if( !params.ref ) error "params.ref (Reference FASTA) must be specified"
if( !params.snp_site_vcf ) error "params.snp_site_vcf (SNP sites VCF) must be specified"
if( !params.indel_site_vcf ) error "params.indel_site_vcf (INDEL sites VCF) must be specified"
if( !params.sv_sites_vcf ) error "params.sv_sites_vcf (SV sites BCF) must be specified"
if( !params.snp_markers_intervals ) error "params.snp_markers_intervals (SNP intervals file) must be specified"
if( !params.indel_markers_intervals ) error "params.indel_markers_intervals (INDEL intervals file) must be specified"
if( !params.gatk ) error "params.gatk_path (GATK jar path) must be specified"
if( !params.picard ) error "params.picard_path (Picard jar path) must be specified"
if( !params.beagle ) error "params.beagle_path (Beagle jar path) must be specified"

workflow {
    bam_ch = Channel.fromPath(params.bam_glob, checkIfExists: true)
    bam_list_ch = bam_ch.collect()

    ref = file(params.ref, checkIfExists: true)
    picard = file(params.picard, checkIfExists: true)
    snp_markers_intervals = file(params.snp_markers_intervals, checkIfExists: true)
    indel_markers_intervals = file(params.indel_markers_intervals, checkIfExists: true)
    snp_site_vcf = file(params.snp_site_vcf, checkIfExists: true)
    indel_site_vcf = file(params.indel_site_vcf, checkIfExists: true)
    sv_sites_vcf = file(params.sv_sites_vcf, checkIfExists: true)

    bam_tuples_ch = bam_ch.map { bam ->
        def sample_id = bam.simpleName
        def bai = file("${bam}.bai", checkIfExists: true)
        tuple(sample_id, bam, bai)
    }

    index_ref_ch = INDEX_REFERENCE(ref, picard)
    snp_vcf_ch = UNIFIED_GENOTYPER_SNP(ref, bam_list_ch, index_ref_ch.fai, index_ref_ch.dict, snp_markers_intervals, snp_site_vcf)
    indel_vcf_ch = UNIFIED_GENOTYPER_INDEL(ref, bam_list_ch, index_ref_ch.fai, index_ref_ch.dict, indel_markers_intervals, indel_site_vcf)

    snp_format_ch = GATK_SNP_FORMAT(snp_vcf_ch.vcf, snp_vcf_ch.tbi, snp_site_vcf)
    indel_format_ch = GATK_INDEL_FORMAT(indel_vcf_ch.vcf_gz, indel_vcf_ch.vcf_gz_index, indel_site_vcf)

    sv_gt_ch = DELLY_SV_GENOTYPE(ref, sv_sites_vcf, bam_tuples_ch)


    bcf_list_ch = sv_gt_ch.map { sample_id, bcf, bcf_index -> bcf}.collect()
    bcf_index_list_ch = sv_gt_ch.map { sample_id, bcf, bcf_index -> bcf_index}.collect()
    
    sv_merged_ch = BCFTOOLS_MERGE_GENOTYPE(bcf_list_ch, bcf_index_list_ch, indel_format_ch.samples_order)

    concat_vcf_ch = CONCAT_VCF(snp_format_ch.vcf, snp_format_ch.tbi, indel_format_ch.vcf_gz, indel_format_ch.vcf_gz_index, sv_merged_ch.vcf_gz, sv_merged_ch.vcf_gz_index)

    beagle_impute_ch = BEAGLE_IMPUTATION(concat_vcf_ch)

    pop_snp_ch = POP_SNP(beagle_impute_ch.impute_vcf)
    pop_indel_ch = POP_INDEL(beagle_impute_ch.impute_vcf)
    pop_sv_ch = POP_SV(beagle_impute_ch.impute_vcf)
}
