process ORTHOFINDER {
    label 'process_high'
    label 'process_long'
    conda "${projectDir}/env/pangenome.yml"
    publishDir "${params.outdir}/06_pangenome", mode: params.publish_dir_mode

    input:
    path faa_dir

    output:
    path "Orthogroups.tsv"                  , emit: orthogroups
    path "Orthogroups.GeneCount.tsv"        , emit: gene_counts
    path "gene_presence_absence.csv"        , emit: presence_absence
    path "single_copy_orthogroups.txt"      , emit: single_copy
    path "pangenome_partitions.tsv"         , emit: partitions
    path "orthofinder_raw"                  , emit: raw
    path "versions.yml"                     , emit: versions

    script:
    """
    orthofinder -f ${faa_dir} -S diamond -og -t ${task.cpus} -o orthofinder_raw

    RES=\$(find orthofinder_raw -type d -name 'Results_*' | head -1)
    cp \$RES/Orthogroups/Orthogroups.tsv .
    cp \$RES/Orthogroups/Orthogroups.GeneCount.tsv .

    partition_pangenome.py \\
        --gene-counts Orthogroups.GeneCount.tsv \\
        --core ${params.core_threshold} \\
        --soft-core ${params.soft_core_threshold} \\
        --shell ${params.shell_threshold} \\
        --out-matrix gene_presence_absence.csv \\
        --out-single-copy single_copy_orthogroups.txt \\
        --out-partitions pangenome_partitions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder -h 2>&1 | grep -m1 -o 'version [0-9.]*' | sed 's/version //')
        diamond: \$(diamond --version 2>&1 | sed 's/.*version //')
    END_VERSIONS
    """
}

process HEAPS_LAW {
    label 'process_low'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/06_pangenome", mode: params.publish_dir_mode

    input:
    path presence_absence

    output:
    path "accumulation_curve.tsv", emit: curve
    path "heaps_law_fit.tsv"     , emit: fit

    script:
    """
    heaps_law.py --matrix ${presence_absence} --permutations 100 \\
        --out-curve accumulation_curve.tsv --out-fit heaps_law_fit.tsv
    """
}
