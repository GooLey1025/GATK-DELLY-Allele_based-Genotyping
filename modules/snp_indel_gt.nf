process INDEX_REFERENCE {
    cpus "${params.threads}"
    input:
    path ref
    path picard

    output:
    path "${ref.name}.fai", emit: fai
    path "${ref.baseName}.dict", emit: dict

    script:
    """
    ${params.samtools} faidx ${ref}
    ${params.java} -jar ${picard} CreateSequenceDictionary R=${ref} O=${ref.baseName}.dict
    """
}

process UNIFIED_GENOTYPER_SNP {
    publishDir "${params.out_dir}/snp_indel_genotype", mode: 'copy'
    cpus "${params.threads}"
    input:
    path ref
    path bam_files
    path fai
    path dict
    path intervals
    path snp_site_vcf
    output:
    path "${params.project}.snp.raw.vcf.gz", emit: vcf
    path "${params.project}.snp.raw.vcf.gz.tbi", emit: tbi

    script:
    def bam_inputs = bam_files.collect { "-I ${it}" }.join(' ')
    """    
    ${params.gatk_java} -Xmx${params.gatk_memory} -jar ${params.gatk} \\
        -R ${ref} \\
        -T UnifiedGenotyper \\
        ${bam_inputs} \\
        -glm SNP \\
        -L ${intervals} \\
        -alleles ${snp_site_vcf} \\
        -gt_mode GENOTYPE_GIVEN_ALLELES \\
        --output_mode EMIT_ALL_SITES \\
        -stand_call_conf 30 \\
        -nct ${params.threads} \\
        -o ${params.project}.snp.raw.vcf.gz

    ${params.tabix} -f -p vcf ${params.project}.snp.raw.vcf.gz
    """
}

process GATK_SNP_FORMAT {
    publishDir "${params.out_dir}/snp_indel_genotype", mode: 'copy'
    cpus params.threads

    input:
    path snp_raw_vcf
    path snp_raw_vcf_index
    path snp_site_vcf

    output:
    path "${params.project}.snp.format.id.vcf.gz", emit: vcf
    path "${params.project}.snp.format.id.vcf.gz.tbi", emit: tbi

    script:
    """
    set -euo pipefail

    tmp_vcf="${params.project}.snp.altfilled.vcf"

    # header
    ${params.bcftools} view -h "${snp_raw_vcf}" > "\$tmp_vcf"

    # fill ALT='.' using sites VCF (match by CHROM,POS,REF)
    awk 'BEGIN{FS=OFS="\\t"}
      NR==FNR{
        key=\$1"\\t"\$2"\\t"\$4
        alt[key]=\$5
        next
      }
      {
        key=\$1"\\t"\$2"\\t"\$4
        if(\$5=="." && (key in alt)) \$5=alt[key]
        print
      }' <(${params.bcftools} view -H "${snp_site_vcf}") <(${params.bcftools} view -H "${snp_raw_vcf}") >> "\$tmp_vcf"

    bash ${params.assign_id} SNP \$tmp_vcf > ${params.project}.snp.format.id.vcf

    # compress + index
    bgzip -@ ${task.cpus} -c "${params.project}.snp.format.id.vcf" > "${params.project}.snp.format.id.vcf.gz"
    ${params.bcftools} index -t --threads ${task.cpus} "${params.project}.snp.format.id.vcf.gz"
    """
}

process UNIFIED_GENOTYPER_INDEL {
    cpus "${params.threads}"
    publishDir "${params.out_dir}/snp_indel_genotype", mode: 'copy'
    input:
    path ref
    path bam_files
    path fai
    path dict
    path intervals
    path indel_site_vcf
    output:
    path "${params.project}.indel.raw.vcf.gz", emit: vcf_gz
    path "${params.project}.indel.raw.vcf.gz.tbi", emit: vcf_gz_index

    script:
    def bam_inputs = bam_files.collect { "-I ${it}" }.join(' ')
    """    
    ${params.gatk_java} -Xmx${params.gatk_memory} -jar ${params.gatk} \\
        -R ${ref} \\
        -T UnifiedGenotyper \\
        ${bam_inputs} \\
        -glm INDEL \\
        -L ${intervals} \\
        -alleles ${indel_site_vcf} \\
        -gt_mode GENOTYPE_GIVEN_ALLELES \\
        --output_mode EMIT_ALL_SITES \\
        -stand_call_conf 30 \\
        -nct ${params.threads} \\
        -o ${params.project}.indel.raw.vcf.gz

    ${params.tabix} -f -p vcf ${params.project}.indel.raw.vcf.gz
    """
}

process GATK_INDEL_FORMAT {
    publishDir "${params.out_dir}/snp_indel_genotype", mode: 'copy'
    cpus params.threads

    input:
    path indel_raw_vcf
    path indel_raw_vcf_index
    path indel_site_vcf

    output:
    path "${params.project}.indel.format.id.vcf.gz", emit: vcf_gz
    path "${params.project}.indel.format.id.vcf.gz.tbi", emit: vcf_gz_index
    path "samples.order.txt", emit: samples_order
    script:
    """
    set -euo pipefail

    tmp_vcf="${params.project}.indel.altfilled.vcf"

    # header
    ${params.bcftools} view -h "${indel_raw_vcf}" > "\$tmp_vcf"

    # fill ALT='.' using sites VCF (match by CHROM,POS,REF)
    awk 'BEGIN{FS=OFS="\\t"}
      NR==FNR{
        key=\$1"\\t"\$2"\\t"\$4
        alt[key]=\$5
        next
      }
      {
        key=\$1"\\t"\$2"\\t"\$4
        if(\$5=="." && (key in alt)) \$5=alt[key]
        print
      }' <(${params.bcftools} view -H "${indel_site_vcf}") <(${params.bcftools} view -H "${indel_raw_vcf}") >> "\$tmp_vcf"

    bash ${params.assign_id} INDEL \$tmp_vcf > ${params.project}.indel.format.id.vcf

    # compress + index
    bgzip -@ ${task.cpus} -c "${params.project}.indel.format.id.vcf" > "${params.project}.indel.format.id.vcf.gz"
    ${params.bcftools} index -t --threads ${task.cpus} "${params.project}.indel.format.id.vcf.gz"
    ${params.bcftools} query -l ${params.project}.indel.format.id.vcf.gz > samples.order.txt
    """
}


