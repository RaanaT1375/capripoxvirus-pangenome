process PARSNP {
    label 'process_high'
    conda "${projectDir}/env/recombination.yml"
    publishDir "${params.outdir}/09_recombination", mode: params.publish_dir_mode,
               pattern: "core_alignment.fasta"

    input:
    path genomes_dir
    path reference

    output:
    path "core_alignment.fasta", emit: alignment
    path "parsnp_raw"          , emit: raw
    path "versions.yml"        , emit: versions

    script:
    // -c forces every input genome into the alignment, so the core block is
    // constrained across all isolates rather than the best-covered subset.
    """
    parsnp -r ${reference} -d ${genomes_dir} -o parsnp_raw -p ${task.cpus} -c
    harvesttools -i parsnp_raw/parsnp.ggr -M core_alignment.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parsnp: \$(parsnp --version 2>&1 | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | head -1)
    END_VERSIONS
    """
}

process GUBBINS {
    label 'process_high'
    label 'process_long'
    conda "${projectDir}/env/recombination.yml"
    publishDir "${params.outdir}/09_recombination", mode: params.publish_dir_mode

    input:
    path core_alignment

    output:
    path "*.recombination_predictions.gff", emit: recombination
    path "*.node_labelled.final_tree.tre" , emit: tree
    path "*.filtered_polymorphic_sites.fasta", emit: snps
    path "rm_ratio.tsv"                   , emit: rm_ratio
    path "versions.yml"                   , emit: versions

    script:
    """
    run_gubbins.py --prefix capripox \\
        --tree-builder ${params.gubbins_treebuilder} \\
        --iterations ${params.gubbins_iterations} \\
        --threads ${task.cpus} \\
        ${core_alignment}

    rm_ratio.py --gff capripox.recombination_predictions.gff \\
        --embl capripox.branch_base_reconstruction.embl \\
        --alignment ${core_alignment} \\
        --out rm_ratio.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: \$(run_gubbins.py --version 2>&1)
    END_VERSIONS
    """
}

process FASTBAPS {
    label 'process_medium'
    conda "${projectDir}/env/recombination.yml"
    publishDir "${params.outdir}/10_population_structure", mode: params.publish_dir_mode

    input:
    path snp_alignment

    output:
    path "fastbaps_clusters.csv", emit: clusters
    path "fastbaps_priors.tsv"  , emit: sensitivity
    path "versions.yml"         , emit: versions

    script:
    // The number of clusters is not fixed a priori; it follows from the
    // Bayesian hierarchical partitioning of the SNP matrix itself.
    """
    run_fastbaps.R ${snp_alignment} fastbaps_clusters.csv fastbaps_priors.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastbaps: \$(Rscript -e 'cat(as.character(packageVersion("fastbaps")))')
    END_VERSIONS
    """
}
