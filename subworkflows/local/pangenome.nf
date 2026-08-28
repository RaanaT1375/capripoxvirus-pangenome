//
// Stages 05-06: proteome-level quality control and pangenome reconstruction.
//

include { DIAMOND_CONTAMINATION } from '../../modules/local/pangenome_qc'
include { PANGENOME_QC_FILTER   } from '../../modules/local/pangenome_qc'
include { ORTHOFINDER           } from '../../modules/local/orthofinder'
include { HEAPS_LAW             } from '../../modules/local/orthofinder'

workflow PANGENOME {

    take:
    ch_assemblies   // [ meta, faa ]
    ch_proteome

    main:
    ch_versions = Channel.empty()

    DIAMOND_CONTAMINATION( ch_assemblies, ch_proteome )
    ch_versions = ch_versions.mix(DIAMOND_CONTAMINATION.out.versions.first())

    PANGENOME_QC_FILTER(
        DIAMOND_CONTAMINATION.out.report.collect(),
        ch_assemblies.map { meta, faa -> faa }.collect(),
        Channel.empty().collect().ifEmpty([])
    )

    ch_retained = PANGENOME_QC_FILTER.out.retained
        .splitText()
        .map { it.trim() }
        .filter { it }

    ch_clean = ch_assemblies
        .map { meta, faa -> [ meta.id, meta, faa ] }
        .join( ch_retained.map { id -> [ id, true ] } )
        .map { id, meta, faa, keep -> [ meta, faa ] }

    ORTHOFINDER( ch_clean.map { meta, faa -> faa }.collect() )
    HEAPS_LAW( ORTHOFINDER.out.presence_absence )
    ch_versions = ch_versions.mix(ORTHOFINDER.out.versions)

    ch_species = ch_clean
        .map { meta, faa -> "${meta.id},${meta.species},${meta.host},${meta.country},${meta.year},${meta.condition}" }
        .collectFile(name: 'species_assignment.csv', newLine: true,
                     seed: 'Name,species,host,country,year,condition',
                     storeDir: "${params.outdir}/05_pangenome_qc")

    emit:
    presence_absence         = ORTHOFINDER.out.presence_absence
    single_copy_orthogroups  = ORTHOFINDER.out.single_copy
    partitions               = ORTHOFINDER.out.partitions
    faa_clean                = ch_clean
    fna_clean                = ch_clean
    species_assignment       = ch_species
    versions                 = ch_versions
}
