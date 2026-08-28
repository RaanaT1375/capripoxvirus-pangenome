process SPADES {
    tag   "${meta.id} (${mode})"
    label 'process_high'
    conda "${projectDir}/env/assembly.yml"
    publishDir "${params.outdir}/03_assembly/${meta.id}", mode: params.publish_dir_mode,
               pattern: "*.{fasta,log}"

    input:
    tuple val(meta), path(reads)
    val   mode          // 'careful' for the first pass, 'rescue' after bacterial depletion

    output:
    tuple val(meta), path("${meta.id}_contigs.fasta"), emit: contigs
    path  "*.log"       , emit: log
    path  "versions.yml", emit: versions

    script:
    // First pass uses --careful. Datasets re-assembled after bacterial
    // depletion use --only-assembler, which bypasses read error correction;
    // this matches the two-pass procedure used in the published analysis.
    def spades_mode = (mode == 'rescue') ? '--only-assembler' : '--careful'
    def minlen      = (mode == 'rescue') ? params.min_contig_len_rescue : params.min_contig_len
    def input_args  = meta.single_end ? "--s1 ${reads[0]}" : "--pe1-1 ${reads[0]} --pe1-2 ${reads[1]}"
    def mem_gb      = task.memory.toGiga()
    """
    spades.py ${input_args} ${spades_mode} \\
        --threads ${task.cpus} --memory ${mem_gb} \\
        -o spades_out > ${meta.id}.spades.log 2>&1

    awk -v L=${minlen} 'BEGIN{RS=">";FS="\\n"} NR>1 {
        seq=""; for(i=2;i<=NF;i++) seq=seq \$i;
        if (length(seq) >= L) printf ">%s\\n%s\\n", \$1, seq
    }' spades_out/contigs.fasta > ${meta.id}_contigs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/SPAdes genome assembler v//')
        min_contig_length: ${minlen}
        mode: ${spades_mode}
    END_VERSIONS
    """
}
