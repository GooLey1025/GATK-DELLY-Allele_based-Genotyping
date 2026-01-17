#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.project = 'output'
params.bam_glob = '/data3/home/gulei/projects/GraphPan/gatk/705rice/bam_dir/chrfix_bam/*.SMfix.bam'
params.out_dir = 'out_dir'


params.ref = 'Nip.chrnum.sorted.fa'
params.snp_site_vcf = '/data3/home/gulei/projects/GraphPan/gwas/markers_d/705rice/705rice/705rice.snp.sites.vcf'
params.indel_site_vcf = '/data3/home/gulei/projects/GraphPan/gwas/markers_d/705rice/705rice/705rice.indel.sites.vcf'
params.snp_markers_intervals = '/data3/home/gulei/projects/GraphPan/gwas/markers_d/705rice/705rice/705rice.snp.markers.intervals'
params.indel_markers_intervals = '/data3/home/gulei/projects/GraphPan/gwas/markers_d/705rice/705rice/705rice.indel.markers.intervals'
params.sv_sites_vcf = '/data3/home/gulei/projects/GraphPan/gwas/markers_d/705rice/705rice/705rice.delly.sv.sites.bcf'

params.threads = 48
params.beagle_path = '/data/home/gulei/softwares/beagle.29Oct24.c8e.jar'
params.gatk_memory = '100g'
params.java_path = '/usr/bin/java' ?: 'java'
params.gatk_path = './GenomeAnalysisTK3.7.jar'
params.picard_path = '/data/home/gulei/softwares/picard.jar' 
params.delly_path = '/data/home/gulei/softwares/delly' ?: 'delly'
params.bcftools_path = 'bcftools'
params.assign_id_path = './scripts/assign_id.sh'

params.assign_id = file(params.assign_id_path).toAbsolutePath() // convert relative path to absolute path
params.gatk = file(params.gatk_path).toAbsolutePath()
params.beagle = file(params.beagle_path).toAbsolutePath()
params.picard = file(params.picard_path).toAbsolutePath()
params.delly = file(params.delly_path).toAbsolutePath()
params.bcftools = file(params.bcftools_path).toAbsolutePath()


params.beagle_memory = '300 GB'
params.beagle_cpus = 64

include { INDEX_REFERENCE; UNIFIED_GENOTYPER_SNP; UNIFIED_GENOTYPER_INDEL; GATK_SNP_FORMAT; GATK_INDEL_FORMAT } from './modules/snp_indel_gt'
include { DELLY_SV_GENOTYPE; BCFTOOLS_MERGE_GENOTYPE } from './modules/sv_gt'

if( !params.gatk ) error "params.gatk (GATK jar path) must be specified"
if( !params.picard ) error "params.picard (Picard jar path) must be specified"

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
    
    sv_merged_ch = BCFTOOLS_MERGE_GENOTYPE(bcf_list_ch, bcf_index_list_ch)

    concat_vcf_ch = CONCAT_VCF(snp_format_ch.vcf, snp_format_ch.tbi, indel_format_ch.vcf_gz, indel_format_ch.vcf_gz_index, sv_merged_ch.vcf_gz, sv_merged_ch.vcf_gz_index)

    beagle_impute_ch = BEAGLE_IMPUTATION(concat_vcf_ch)

    pop_snp_ch = POP_SNP(beagle_impute_ch.impute_vcf)
    pop_indel_ch = POP_INDEL(beagle_impute_ch.impute_vcf)
    pop_sv_ch = POP_SV(beagle_impute_ch.impute_vcf)
}


process CONCAT_VCF {
    input:
    path snp_vcf_gz
    path snp_vcf_gz_index
    path indel_vcf_gz
    path indel_vcf_gz_index
    path sv_vcf_gz
    path sv_vcf_gz_index
    output:
    path "${params.project}.snp.indel.sv.vcf"

    script:
    """
    bcftools concat -a ${snp_vcf_gz} ${indel_vcf_gz} ${sv_vcf_gz} -o ${params.project}.snp.indel.sv.vcf
    """
}

process BEAGLE_IMPUTATION {
    publishDir params.out_dir, mode: 'copy', pattern: ".vcf.gz"
    memory "${params.beagle_memory}"
    cpus "${params.beagle_cpus}"
    input:
    path snp_indel_sv_vcf

    output:
    path "${params.project}.snp.indel.sv.impute.vcf.gz", emit: impute_vcf

    script:
    """
    mkdir -p ./beagle_TMP
    java -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./beagle_TMP \\
        -jar ${params.beagle} \\
        gt=${snp_indel_sv_vcf} \\
        out=./${params.project}.snp.indel.sv.impute
    """
}

process POP_SNP {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_indel_sv_impute_vcf
    output:
    path "${params.project}.snp.impute.vcf.gz", emit: snp_vcf
    script:
    """
    bcftools view -i 'ID ~ "^SNP-"' ${snp_indel_sv_impute_vcf} -o ${params.project}.snp.impute.vcf.gz
    """
}
process POP_INDEL {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_indel_sv_impute_vcf
    output:
    path "${params.project}.indel.impute.vcf.gz", emit: indel_vcf
    script:
    """
    bcftools view -i 'ID ~ "^INDEL-"' ${snp_indel_sv_impute_vcf} -o ${params.project}.indel.impute.vcf.gz
    """
}
process POP_SV {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_indel_sv_impute_vcf
    output:
    path "${params.project}.sv.impute.vcf.gz", emit: sv_vcf
    script:
    """
    bcftools view -i 'ID ~ "^SV-"' ${snp_indel_sv_impute_vcf} -Oz -o ${params.project}.sv.impute.vcf.gz
    """
}