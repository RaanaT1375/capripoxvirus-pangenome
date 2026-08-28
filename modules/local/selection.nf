process PAL2NAL {
    tag   "${og}"
    label 'process_single'
    conda "${projectDir}/env/selection.yml"

    input:
    tuple val(og), path(aa_aln), path(cds)

    output:
    tuple val(og), path("${og}_codon.fasta"), emit: codon
    path  "versions.yml"                    , emit: versions

    script:
    """
    pal2nal.pl ${aa_aln} ${cds} -output fasta -nogap > ${og}_codon.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pal2nal: 14
    END_VERSIONS
    """
}

process HYPHY_BUSTED {
    tag   "${og}"
    label 'hyphy'
    conda "${projectDir}/env/selection.yml"
    publishDir "${params.outdir}/14_selection/busted", mode: params.publish_dir_mode

    input:
    tuple val(og), path(codon_aln), path(tree)

    output:
    tuple val(og), path("${og}_BUSTED.json"), emit: json, optional: true
    path  "${og}_BUSTED.log"                , emit: log
    path  "versions.yml"                    , emit: versions

    script:
    // Synonymous rate variation is active in the HyPhy default configuration.
    // TOLERATE_NUMERICAL_ERRORS is applied on retry only: a handful of highly
    // conserved loci fail optimisation under the default tolerances.
    def tolerate = task.attempt > 1 ? 'ENV=TOLERATE_NUMERICAL_ERRORS=1' : ''
    """
    hyphy busted ${tolerate} \\
        --alignment ${codon_aln} \\
        --tree ${tree} \\
        --output ${og}_BUSTED.json \\
        CPU=${task.cpus} > ${og}_BUSTED.log 2>&1 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hyphy: \$(hyphy --version 2>&1 | grep -o '[0-9.]*' | head -1)
    END_VERSIONS
    """
}

process HYPHY_RELAX {
    tag   "${og}"
    label 'hyphy'
    conda "${projectDir}/env/selection.yml"
    publishDir "${params.outdir}/14_selection/relax", mode: params.publish_dir_mode

    input:
    tuple val(og), path(codon_aln), path(labelled_tree)

    output:
    tuple val(og), path("${og}_RELAX.json"), emit: json, optional: true
    path  "${og}_RELAX.log"                , emit: log

    script:
    """
    hyphy relax --alignment ${codon_aln} --tree ${labelled_tree} \\
        --test Test --reference Reference \\
        --models ${params.relax_models} \\
        --output ${og}_RELAX.json \\
        CPU=${task.cpus} > ${og}_RELAX.log 2>&1 || true
    """
}

process HYPHY_ABSREL {
    tag   "${og}"
    label 'hyphy'
    conda "${projectDir}/env/selection.yml"
    publishDir "${params.outdir}/14_selection/absrel", mode: params.publish_dir_mode

    input:
    tuple val(og), path(codon_aln), path(labelled_tree)

    output:
    tuple val(og), path("${og}_ABSREL.json"), emit: json, optional: true
    path  "${og}_ABSREL.log"                , emit: log

    script:
    """
    hyphy absrel --alignment ${codon_aln} --tree ${labelled_tree} \\
        --branches Stem \\
        --output ${og}_ABSREL.json \\
        CPU=${task.cpus} > ${og}_ABSREL.log 2>&1 || true
    """
}

process SELECTION_DIAGNOSTICS {
    label 'process_low'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/14_selection", mode: params.publish_dir_mode

    input:
    path busted_jsons
    path codon_alignments

    output:
    path "busted_summary.csv"     , emit: summary
    path "dnds_applicability.tsv" , emit: diagnostics

    script:
    // Records, per gene, the number of distinct protein haplotypes, the
    // proportion of sites with dS approximately zero, and whether the
    // likelihood-ratio statistic degenerated to zero. These diagnostics
    // determine whether codon models are interpretable for the dataset at all.
    """
    busted_summary.py --json-dir . --alignments . \\
        --fdr ${params.hyphy_fdr} \\
        --out-summary busted_summary.csv \\
        --out-diagnostics dnds_applicability.tsv
    """
}

process MK_TEST {
    label 'process_medium'
    conda "${projectDir}/env/python.yml"
    publishDir "${params.outdir}/15_mk_test", mode: params.publish_dir_mode

    input:
    path codon_alignments
    path species_assignment
    path traits

    output:
    path "mk_polarized_per_gene.csv"  , emit: per_gene
    path "mk_genome_wide.tsv"         , emit: genome_wide
    path "mk_threshold_sweep.csv"     , emit: sweep
    path "apobec_diagnostic.tsv"      , emit: apobec

    script:
    // Polarised MK test: the ancestral codon is the consensus shared by the
    // SPPV and GTPV outgroups, fixed differences are assigned to the LSDV
    // branch, and codons that are both polymorphic and divergent are counted
    // as polymorphism only.
    """
    # The polarised MK implementation is the archived analysis code, run
    # verbatim so that the published estimate is reproduced exactly. The
    # directory layout it expects is staged here rather than rewritten.
    mkdir -p 10_Selection_Pressure/01_codon_alignments 14_MK_Test/02_results \\
             00_Metadata 11_Scoary/01_inputs
    cp *_codon.fasta 10_Selection_Pressure/01_codon_alignments/
    cp ${species_assignment} 00_Metadata/species_assignment.csv
    cp ${traits} 11_Scoary/01_inputs/traits.csv

    python3 ${projectDir}/bin/legacy/mk_polarized.py
    python3 ${projectDir}/bin/legacy/mk_threshold_sweep.py || true
    python3 ${projectDir}/bin/legacy/apobec_diagnostic.py || true

    cp 14_MK_Test/02_results/mk_polarized_LSDV_branch.csv mk_polarized_per_gene.csv
    cp 14_MK_Test/03_summary/mk_threshold_sweep.csv       mk_threshold_sweep.csv
    cp 14_MK_Test/03_summary/apobec_diagnostic.tsv        apobec_diagnostic.tsv 2>/dev/null || touch apobec_diagnostic.tsv
    grep -m1 alpha 14_MK_Test/02_results/*.csv > mk_genome_wide.tsv || touch mk_genome_wide.tsv
    """
}
