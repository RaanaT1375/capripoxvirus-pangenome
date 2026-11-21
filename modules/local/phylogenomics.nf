process MAFFT_ALIGN {
    tag   "${og}"
    label 'process_low'
    conda "${projectDir}/env/phylogenomics.yml"

    input:
    tuple val(og), path(faa)

    output:
    tuple val(og), path("${og}.aln.faa"), emit: alignment
    path  "versions.yml"                , emit: versions

    script:
    """
    mafft --auto --thread ${task.cpus} ${faa} > ${og}.aln.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//;s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    printf '>t1\\nMSTN\\n>t2\\nMSTN\\n' > ${og}.aln.faa
    touch versions.yml
    """
}

process CONCAT_SUPERMATRIX {
    label 'process_low'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/07_supermatrix", mode: params.publish_dir_mode

    input:
    path alignments

    output:
    path "supermatrix.faa"      , emit: supermatrix
    path "partitions.nex"       , emit: partitions
    path "supermatrix_stats.tsv", emit: stats

    script:
    """
    concat_alignments.py --alignments ${alignments} \\
        --out-fasta supermatrix.faa \\
        --out-partitions partitions.nex \\
        --out-stats supermatrix_stats.tsv
    """

    stub:
    """
    printf '>t1\\nMSTN\\n>t2\\nMSTN\\n' > supermatrix.faa
    printf '#nexus\\nbegin sets;\\nend;\\n' > partitions.nex
    touch supermatrix_stats.tsv
    """
}

process IQTREE {
    label 'process_high'
    label 'process_long'
    conda "${projectDir}/env/phylogenomics.yml"
    publishDir "${params.outdir}/08_phylogeny", mode: params.publish_dir_mode

    input:
    path supermatrix
    path partitions

    output:
    path "core_genome_ML.contree", emit: tree
    path "core_genome_ML.*"      , emit: all
    path "versions.yml"          , emit: versions

    script:
    """
    iqtree -s ${supermatrix} -p ${partitions} \\
        -m ${params.iqtree_model} -bb ${params.iqtree_bootstrap} \\
        -nt AUTO -ntmax ${task.cpus} \\
        -pre core_genome_ML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree --version 2>&1 | grep -m1 -o 'version [0-9.]*' | sed 's/version //')
    END_VERSIONS
    """

    stub:
    """
    printf '(t1:0.1,t2:0.1);\\n' > core_genome_ML.contree
    touch core_genome_ML.log
    touch versions.yml
    """
}
