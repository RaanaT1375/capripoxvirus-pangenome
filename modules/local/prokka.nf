process PROKKA {
    tag   "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/env/annotation.yml"
    publishDir "${params.outdir}/04_annotation/${meta.id}", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(contigs)
    path  reference_proteome

    output:
    tuple val(meta), path("${meta.id}.faa"), emit: faa
    tuple val(meta), path("${meta.id}.ffn"), emit: ffn
    tuple val(meta), path("${meta.id}.gff"), emit: gff
    tuple val(meta), path("${meta.id}.tsv"), emit: tsv
    path  "versions.yml"                   , emit: versions

    script:
    // Homology-directed annotation against the LSDV RefSeq proteome, with the
    // viral kingdom parameter to prevent the over-annotation seen when the
    // bacterial defaults are used on viral contigs.
    """
    prokka --outdir prokka_out --prefix ${meta.id} \\
        --kingdom Viruses \\
        --proteins ${reference_proteome} \\
        --evalue ${params.diamond_evalue} \\
        --cpus ${task.cpus} --force ${contigs} > /dev/null 2>&1

    for ext in faa ffn gff tsv; do cp prokka_out/${meta.id}.\$ext .; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(prokka --version 2>&1 | sed 's/prokka //')
        annotation_mode: with_reference_proteome
    END_VERSIONS
    """
}
