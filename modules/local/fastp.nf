process FASTP {
    tag   "${meta.id}"
    label 'process_low'
    conda "${projectDir}/env/preprocessing.yml"
    publishDir "${params.outdir}/01_read_qc", mode: params.publish_dir_mode,
               pattern: "*.{json,html}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.json")         , emit: json
    path  "*.html"                          , emit: html
    path  "versions.yml"                    , emit: versions

    script:
    def common = """--qualified_quality_phred ${params.fastp_qual_phred} \\
        --unqualified_percent_limit ${params.fastp_unqualified_percent} \\
        --cut_front --cut_tail --cut_mean_quality ${params.fastp_cut_mean_quality} \\
        --length_required ${params.fastp_min_length} \\
        --trim_poly_g --trim_poly_x \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html"""
    if (meta.single_end)
        """
        fastp --in1 ${reads[0]} --out1 ${meta.id}.trim.fastq.gz ${common}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    else
        """
        fastp --in1 ${reads[0]} --in2 ${reads[1]} \\
              --out1 ${meta.id}_1.trim.fastq.gz --out2 ${meta.id}_2.trim.fastq.gz \\
              --detect_adapter_for_pe ${common}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
}
