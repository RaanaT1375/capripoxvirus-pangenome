process BOWTIE2_HOST_DEPLETION {
    tag   "${meta.id} (${meta.host})"
    label 'process_medium'
    conda "${projectDir}/env/preprocessing.yml"
    publishDir "${params.outdir}/02_host_depletion/logs", mode: params.publish_dir_mode,
               pattern: "*.log"

    input:
    tuple val(meta), path(reads)
    path  index_dir

    output:
    tuple val(meta), path("*host_removed*.fastq.gz"), emit: reads
    path  "*.log"       , emit: log
    path  "versions.yml", emit: versions

    script:
    // Host index is selected from the isolate metadata: Bos taurus for cattle,
    // Ovis aries for sheep, Capra hircus for goat.
    def idx = [ 'Bos taurus':'Bos_taurus', 'Ovis aries':'Ovis_aries',
                'Capra hircus':'Capra_hircus' ].get(meta.host, 'Bos_taurus')
    if (meta.single_end)
        """
        bowtie2 -p ${task.cpus} -x ${index_dir}/${idx} \\
                -U ${reads[0]} \\
                --un-gz ${meta.id}_host_removed.fastq.gz \\
                > /dev/null 2> ${meta.id}.bowtie2.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*version //')
        END_VERSIONS
        """
    else
        """
        bowtie2 -p ${task.cpus} -x ${index_dir}/${idx} \\
                -1 ${reads[0]} -2 ${reads[1]} \\
                --un-conc-gz ${meta.id}_host_removed_R%.fastq.gz \\
                > /dev/null 2> ${meta.id}.bowtie2.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*version //')
        END_VERSIONS
        """
}
