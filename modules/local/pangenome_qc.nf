process DIAMOND_CONTAMINATION {
    tag   "${meta.id}"
    label 'process_low'
    conda "${projectDir}/env/annotation.yml"

    input:
    tuple val(meta), path(faa)
    path  reference_proteome

    output:
    path "*.contamination.tsv", emit: report
    path "versions.yml"       , emit: versions

    script:
    """
    diamond makedb --in ${reference_proteome} -d ref --quiet
    diamond blastp -q ${faa} -d ref -o hits.tsv \\
        --evalue ${params.diamond_evalue} --max-target-seqs 1 \\
        --threads ${task.cpus} --quiet

    TOTAL=\$(grep -c '^>' ${faa})
    HIT=\$(cut -f1 hits.tsv | sort -u | wc -l)
    printf "%s\\t%s\\t%s\\t%.4f\\n" "${meta.id}" "\$TOTAL" "\$HIT" \\
        "\$(echo "1 - \$HIT / \$TOTAL" | bc -l)" > ${meta.id}.contamination.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
}

process PANGENOME_QC_FILTER {
    label 'process_single'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/05_pangenome_qc", mode: params.publish_dir_mode

    input:
    path contamination_reports
    path faa_files
    path cds_counts

    output:
    path "retained_genomes.txt" , emit: retained
    path "excluded_genomes.tsv" , emit: excluded
    path "qc_summary.tsv"       , emit: summary

    script:
    """
    pangenome_qc.py \\
        --contamination ${contamination_reports} \\
        --faa-dir . \\
        --cds-counts ${cds_counts} \\
        --contamination-max ${params.contamination_max} \\
        --near-core-threshold ${params.near_core_threshold} \\
        --max-missing-near-core ${params.max_missing_near_core} \\
        --out-retained retained_genomes.txt \\
        --out-excluded excluded_genomes.tsv \\
        --out-summary  qc_summary.tsv
    """
}
