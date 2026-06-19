process CONCAT_VCF {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_vcf_gz
    path snp_vcf_gz_index
    path indel_vcf_gz
    path indel_vcf_gz_index
    path sv_vcf_gz
    path sv_vcf_gz_index
    output:
    path "${params.project}.snp.indel.sv.vcf.gz"

    script:
    """
    ${params.bcftools} concat -a ${snp_vcf_gz} ${indel_vcf_gz} ${sv_vcf_gz} -o ${params.project}.snp.indel.sv.vcf.gz
    """
}

process BEAGLE_IMPUTATION {
    publishDir params.out_dir, mode: 'copy', pattern: ".vcf.gz"
    memory "${params.beagle_memory}"
    cpus "${params.beagle_cpus}"
    input:
    path snp_indel_sv_vcf

    output:
    path "${params.project}.snp.indel.sv.impute.biallelic.vcf.gz", emit: impute_biallelic_vcf

    script:
    """
    mkdir -p ./beagle_TMP
    ${params.java} -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=./beagle_TMP \\
        -jar ${params.beagle} \\
        gt=${snp_indel_sv_vcf} \\
        out=${params.project}.snp.indel.sv.impute

    ${params.bcftools} view -h ${params.project}.snp.indel.sv.impute.vcf.gz | head -n -1 > ${params.project}.snp.indel.sv.header
    ${params.bcftools} view -h ${snp_indel_sv_vcf} >> ${params.project}.snp.indel.sv.header
    ${params.bgzip} ${params.project}.snp.indel.sv.impute.vcf.gz -d -c > ${params.project}.snp.indel.sv.impute.vcf
    ${params.bcftools} reheader -h ${params.project}.snp.indel.sv.header ${params.project}.snp.indel.sv.impute.vcf > ${params.project}.snp.indel.sv.impute.reheader.vcf

    ${params.bcftools} norm -m -both ${params.project}.snp.indel.sv.impute.reheader.vcf -o ${params.project}.snp.indel.sv.impute.biallelic.vcf.gz
    """
}

process POP_SNP {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_indel_sv_impute_biallelic_vcf
    output:
    path "${params.project}.snp.impute.biallelic.vcf.gz", emit: snp_vcf
    script:
    """
    ${params.bcftools} view -i 'ID ~ "^SNP-"' ${snp_indel_sv_impute_biallelic_vcf} -o ${params.project}.snp.impute.biallelic.vcf.gz
    """
}

process POP_INDEL {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_indel_sv_impute_biallelic_vcf
    output:
    path "${params.project}.indel.impute.biallelic.vcf.gz", emit: indel_vcf
    script:
    """
    ${params.bcftools} view -i 'ID ~ "^INDEL-"' ${snp_indel_sv_impute_biallelic_vcf} -o ${params.project}.indel.impute.biallelic.vcf.gz
    """
}

process POP_SV {
    publishDir params.out_dir, mode: 'copy', pattern: "*.vcf.gz"
    input:
    path snp_indel_sv_impute_biallelic_vcf
    output:
    path "${params.project}.sv.impute.biallelic.vcf.gz", emit: sv_vcf
    script:
    """
    ${params.bcftools} view -i 'ID ~ "^SV-"' ${snp_indel_sv_impute_biallelic_vcf} -Oz -o ${params.project}.sv.impute.biallelic.vcf.gz
    """
}
