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
    ${params.bcftools} concat -a ${snp_vcf_gz} ${indel_vcf_gz} ${sv_vcf_gz} -o ${params.project}.snp.indel.sv.vcf
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
    ${params.java} -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./beagle_TMP \\
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
    ${params.bcftools} view -i 'ID ~ "^SNP-"' ${snp_indel_sv_impute_vcf} -o ${params.project}.snp.impute.vcf.gz
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
    ${params.bcftools} view -i 'ID ~ "^INDEL-"' ${snp_indel_sv_impute_vcf} -o ${params.project}.indel.impute.vcf.gz
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
    ${params.bcftools} view -i 'ID ~ "^SV-"' ${snp_indel_sv_impute_vcf} -Oz -o ${params.project}.sv.impute.vcf.gz
    """
}
