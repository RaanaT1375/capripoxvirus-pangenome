process SCOARY {
    label 'process_medium'
    conda "${projectDir}/env/pangwas.yml"
    publishDir "${params.outdir}/13_pangwas/raw", mode: params.publish_dir_mode

    input:
    path presence_absence
    path traits
    path tree

    output:
    path "*.results.csv", emit: results
    path "versions.yml" , emit: versions

    script:
    """
    scoary -g ${presence_absence} -t ${traits} --newicktree ${tree} \\
        -e ${params.scoary_permutations} \\
        --start_col ${params.scoary_start_col} \\
        --threads ${task.cpus} --no-time -o .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scoary: \$(scoary --version 2>&1 | grep -o '[0-9.]*' | head -1)
    END_VERSIONS
    """
}

process SCOARY_SUMMARISE {
    label 'process_single'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/13_pangwas", mode: params.publish_dir_mode

    input:
    path scoary_results
    path traits
    path presence_absence

    output:
    path "S6_pangwas_per_trait_summary.csv", emit: summary
    path "significant_pairs.csv"           , emit: pairs

    script:
    // An association is retained only if it satisfies BOTH the per-trait
    // Bonferroni correction and the permutation empirical p. Tier 1
    // (convergent) associations additionally require >= 3 independent
    // supporting pairs with net positive support.
    """
    scoary_summarise.py --results-dir . --traits ${traits} \\
        --matrix ${presence_absence} \\
        --bonferroni ${params.scoary_bonferroni_p} \\
        --empirical ${params.scoary_empirical_p} \\
        --min-pairs ${params.tier1_min_pairs} \\
        --out-summary S6_pangwas_per_trait_summary.csv \\
        --out-pairs   significant_pairs.csv
    """
}
