process CDS_SCREEN {
    tag   "${meta.id}"
    label 'process_single'
    conda "${projectDir}/env/python.yml"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), env(FLAG), path(tsv), emit: verdict
    path  "*.cds_count.tsv"              , emit: counts

    script:
    // Assemblies carrying substantially more predicted CDS than the ~156 ORFs
    // expected of a capripoxvirus genome are flagged for bacterial depletion
    // and re-assembly, rather than being filtered blindly at the read level.
    """
    N=\$(grep -c \$'\\tCDS\\t' ${tsv} || echo 0)
    if [ "\$N" -gt ${params.cds_flag_threshold} ]; then FLAG=flagged; else FLAG=clean; fi
    printf "%s\\t%s\\t%s\\n" "${meta.id}" "\$N" "\$FLAG" > ${meta.id}.cds_count.tsv
    """

    stub:
    """
    FLAG=clean
    printf '%s\\t1\\tclean\\n' "${meta.id}" > ${meta.id}.cds_count.tsv
    """
}
