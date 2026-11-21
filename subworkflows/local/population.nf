//
// Stages 09-13: recombination, population structure, phylogeography,
// temporal signal and pan-GWAS.
//

include { PARSNP           } from '../../modules/local/recombination'
include { GUBBINS          } from '../../modules/local/recombination'
include { FASTBAPS         } from '../../modules/local/recombination'
include { PASTML           } from '../../modules/local/phylogeography'
include { SUBSAMPLE_TREES  } from '../../modules/local/phylogeography'
include { ROOT_TO_TIP      } from '../../modules/local/phylogeography'
include { SCOARY           } from '../../modules/local/scoary'
include { SCOARY_SUMMARISE } from '../../modules/local/scoary'

workflow POPULATION {

    take:
    ch_genomes          // [ meta, fasta ]
    ch_tree             // core-genome ML tree
    ch_presence_absence
    ch_reference

    main:
    ch_versions = Channel.empty()

    // ---- 09 recombination --------------------------------------------------
    PARSNP( ch_genomes.map { meta, f -> f }.collect(), ch_reference )
    GUBBINS( PARSNP.out.alignment )
    ch_versions = ch_versions.mix(PARSNP.out.versions, GUBBINS.out.versions)

    // ---- 10 population structure ------------------------------------------
    FASTBAPS( GUBBINS.out.snps )
    ch_versions = ch_versions.mix(FASTBAPS.out.versions)

    // ---- metadata assembled from the samplesheet ---------------------------
    ch_metadata = ch_genomes
        .map { meta, f -> "${meta.id}\t${meta.country}\t${meta.year}\t${meta.species}\t${meta.host}\t${meta.condition}" }
        .collectFile(name: 'isolate_metadata.tsv', newLine: true,
                     seed: 'id\tContinent\tYear\tSpecies\tHost\tCondition',
                     storeDir: "${params.outdir}/11_phylogeography")

    // ---- 11 phylogeography, with a balanced sub-sampling control -----------
    if (!params.skip_phylogeography) {
        PASTML( Channel.of('full').combine(GUBBINS.out.tree).combine(ch_metadata) )
        SUBSAMPLE_TREES( GUBBINS.out.tree, ch_metadata )

        ch_reps = SUBSAMPLE_TREES.out.trees
            .flatten()
            .map { t -> [ t.baseName.replaceAll('_tree$',''), t ] }
            .join( SUBSAMPLE_TREES.out.metadata.flatten()
                     .map { m -> [ m.baseName.replaceAll('_metadata$',''), m ] } )

        PASTML( ch_reps )
        ch_versions = ch_versions.mix(PASTML.out.versions.first())
    }

    // ---- 12 temporal signal ------------------------------------------------
    ROOT_TO_TIP( GUBBINS.out.tree, ch_metadata )

    // ---- 13 pan-GWAS -------------------------------------------------------
    ch_traits = ch_genomes
        .map { meta, f -> meta }
        .collectFile(name: 'traits.csv') { meta ->
            // Binary trait matrix built from the samplesheet metadata.
            def africa = meta.country == 'Africa' ? 1 : 0
            def asia   = meta.country == 'Asia'   ? 1 : 0
            def europe = meta.country == 'Europe' ? 1 : 0
            def vacc   = meta.condition == 'vaccine' ? 1 : 0
            "${meta.id},${africa},${asia},${europe},${vacc}\n"
        }

    SCOARY( ch_presence_absence, ch_traits, ch_tree )
    SCOARY_SUMMARISE( SCOARY.out.results.collect(), ch_traits, ch_presence_absence )
    ch_versions = ch_versions.mix(SCOARY.out.versions)

    emit:
    recombination = GUBBINS.out.recombination
    rm_ratio      = GUBBINS.out.rm_ratio
    clusters      = FASTBAPS.out.clusters
    gwas_summary  = SCOARY_SUMMARISE.out.summary
    metadata      = ch_metadata
    versions      = ch_versions
}
