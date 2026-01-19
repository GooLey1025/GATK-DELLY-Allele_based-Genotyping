process DELLY_SV_GENOTYPE {

    tag "${sample_id}"
    publishDir "${params.out_dir}/sample_genotype", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    path ref
    path sv_sites_vcf
    tuple val(sample_id), path(bam_file), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}.geno.bcf"), path("${sample_id}.geno.bcf.csi")

    script:
    """
    ${params.delly} call -g ${ref} -v ${sv_sites_vcf} -o ${sample_id}.geno.bcf ${bam_file} > ${sample_id}.geno.log 2>&1

    nvar=\$(${params.bcftools} view -H ${sample_id}.geno.bcf 2>/dev/null | wc -l | awk '{print \$1}')
    if [ "\$nvar" -eq 0 ]; then
        echo "[FATAL] DELLY produced 0 variant records for sample ${sample_id}" >&2
        echo "[FATAL] Check input sites VCF and BAM mapping/reference consistency." >&2
        exit 1
    fi

    ${params.bcftools} index -f ${sample_id}.geno.bcf
    """
}

process BCFTOOLS_MERGE_GENOTYPE {

    cpus 8
    memory '32 GB'
    publishDir "${params.out_dir}/sv_genotype", mode: 'copy'
    
    input:
    path genotype_bcfs
    path genotype_bcfs_index
    path samples_order
    output:
    path "${params.project}.sv.id.SampleSort.vcf.gz", emit: vcf_gz
    path "${params.project}.sv.id.SampleSort.vcf.gz.csi", emit: vcf_gz_index

    script:
    """
    ${params.bcftools} merge --threads ${task.cpus} -m id -o ${params.project}.sv.vcf ${genotype_bcfs.join(' ')}
    bash ${params.assign_id} SV ${params.project}.sv.vcf > ${params.project}.sv.id.vcf
    ${params.bcftools} view --threads ${task.cpus} -S ${samples_order} ${params.project}.sv.id.vcf -W -o ${params.project}.sv.id.SampleSort.vcf.gz
    """
}

