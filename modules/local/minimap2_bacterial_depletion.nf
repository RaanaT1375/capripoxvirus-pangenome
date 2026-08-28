process MINIMAP2_BACTERIAL_DEPLETION {
    tag   "${meta.id}"
    conda "${projectDir}/env/preprocessing.yml"
    publishDir "${params.outdir}/02_bacterial_depletion", mode: params.publish_dir_mode,
               pattern: "*.stats"

    input:
    tuple val(meta), path(reads)
    path  bacteria_mmi

    output:
    tuple val(meta), path("*pure_viral*.fastq.gz"), emit: reads
    path  "*.stats"     , emit: stats
    path  "versions.yml", emit: versions

    script:
    // Negative selection: a read pair is retained only when BOTH mates fail to
    // align to the bacterial reference set (samtools view -f 12). Mapping the
    // reads positively against the LSDV reference was deliberately avoided,
    // since that biases recovered gene content towards the reference and
    // depletes the divergent terminal regions.
    def args = task.ext.args ?: '-ax sr'
    if (meta.single_end)
        """
        minimap2 ${args} -t ${task.cpus} ${bacteria_mmi} ${reads[0]} \\
            | samtools view -b -f 4 - \\
            | samtools fastq - | gzip > ${meta.id}_pure_viral.fastq.gz

        echo "sample\tretained_reads" > ${meta.id}.depletion.stats
        echo -e "${meta.id}\t\$(zcat ${meta.id}_pure_viral.fastq.gz | wc -l | awk '{print \$1/4}')" \\
            >> ${meta.id}.depletion.stats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version)
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
    else
        """
        minimap2 ${args} -t ${task.cpus} ${bacteria_mmi} ${reads[0]} ${reads[1]} \\
            | samtools view -b -f 12 - \\
            | samtools sort -n -@ ${task.cpus} - \\
            | samtools fastq -1 ${meta.id}_pure_viral_R1.fastq -2 ${meta.id}_pure_viral_R2.fastq -0 /dev/null -s /dev/null -
        gzip ${meta.id}_pure_viral_R1.fastq ${meta.id}_pure_viral_R2.fastq

        echo -e "sample\tretained_pairs" > ${meta.id}.depletion.stats
        echo -e "${meta.id}\t\$(zcat ${meta.id}_pure_viral_R1.fastq.gz | wc -l | awk '{print \$1/4}')" \\
            >> ${meta.id}.depletion.stats

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version)
            samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        END_VERSIONS
        """
}
