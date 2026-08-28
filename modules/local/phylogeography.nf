process PASTML {
    tag   "${label}"
    label 'process_medium'
    conda "${projectDir}/env/phylogeography.yml"
    publishDir "${params.outdir}/11_phylogeography/${label}", mode: params.publish_dir_mode

    input:
    tuple val(label), path(tree), path(metadata)

    output:
    tuple val(label), path("*/marginal_probabilities.character_Continent*.tab"), emit: marginal
    tuple val(label), path("*/combined_ancestral_states.tab")                  , emit: states
    path  "versions.yml"                                                        , emit: versions

    script:
    """
    pastml --tree ${tree} --data ${metadata} \\
        --columns Continent \\
        --model ${params.pastml_model} \\
        --prediction_method ${params.pastml_prediction} \\
        --work_dir ${label}_pastml \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pastml: \$(pastml --version 2>&1 | grep -o '[0-9.]*' | head -1)
    END_VERSIONS
    """
}

process SUBSAMPLE_TREES {
    label 'process_single'
    conda "${projectDir}/env/python.yml"

    input:
    path tree
    path metadata

    output:
    path "rep*_tree.nwk"    , emit: trees
    path "rep*_metadata.tsv", emit: metadata

    script:
    // Continental representation is equalised to the least-sampled continent
    // so that the inferred root state can be tested against sampling bias.
    """
    balanced_subsample.py --tree ${tree} --metadata ${metadata} \\
        --replicates ${params.subsampling_replicates} \\
        --column Continent --seed ${params.mk_seed}
    """
}

process ROOT_TO_TIP {
    label 'process_single'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/12_temporal_signal", mode: params.publish_dir_mode

    input:
    path tree
    path metadata

    output:
    path "root_to_tip_distances.tsv", emit: distances
    path "regression_summary.tsv"   , emit: regression

    script:
    """
    root_to_tip.py --tree ${tree} --metadata ${metadata} \\
        --out-distances root_to_tip_distances.tsv \\
        --out-regression regression_summary.tsv
    """
}
