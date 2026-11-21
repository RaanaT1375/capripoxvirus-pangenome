//
// Stages 07-08: per-orthogroup alignment, supermatrix, partitioned ML tree.
//

include { MAFFT_ALIGN        } from '../../modules/local/phylogenomics'
include { CONCAT_SUPERMATRIX } from '../../modules/local/phylogenomics'
include { IQTREE             } from '../../modules/local/phylogenomics'

workflow PHYLOGENOMICS {

    take:
    ch_single_copy   // text file listing single-copy orthogroup IDs
    ch_faa           // [ meta, faa ]

    main:
    ch_versions = Channel.empty()

    ch_og = ch_single_copy
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { og -> [ og, file("${params.outdir}/06_pangenome/orthogroup_sequences/${og}.fa") ] }

    MAFFT_ALIGN( ch_og )
    ch_versions = ch_versions.mix(MAFFT_ALIGN.out.versions.first())

    CONCAT_SUPERMATRIX( MAFFT_ALIGN.out.alignment.map { og, aln -> aln }.collect() )

    IQTREE( CONCAT_SUPERMATRIX.out.supermatrix, CONCAT_SUPERMATRIX.out.partitions )
    ch_versions = ch_versions.mix(IQTREE.out.versions)

    emit:
    tree            = IQTREE.out.tree
    supermatrix     = CONCAT_SUPERMATRIX.out.supermatrix
    gene_alignments = MAFFT_ALIGN.out.alignment
    versions        = ch_versions
}
