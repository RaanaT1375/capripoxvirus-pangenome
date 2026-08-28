//
// Stages 14-16: codon-model selection scans, their diagnostics, and the
// polarised McDonald-Kreitman test.
//
// The codon-model scans are run for two reasons: to test the genus-level
// hypothesis, and to establish whether dN/dS models are applicable within the
// clonal LSDV population at all. SELECTION_DIAGNOSTICS records the evidence for
// the latter, which is what motivates the population-genetic framework.
//

include { PAL2NAL               } from '../../modules/local/selection'
include { HYPHY_BUSTED          } from '../../modules/local/selection'
include { HYPHY_RELAX           } from '../../modules/local/selection'
include { HYPHY_ABSREL          } from '../../modules/local/selection'
include { SELECTION_DIAGNOSTICS } from '../../modules/local/selection'
include { MK_TEST               } from '../../modules/local/selection'

workflow SELECTION {

    take:
    ch_alignments        // [ og, aa_alignment ]
    ch_tree
    ch_species

    main:
    ch_versions = Channel.empty()

    ch_cds = ch_alignments.map { og, aln ->
        [ og, aln, file("${params.outdir}/06_pangenome/orthogroup_cds/${og}.ffn") ]
    }

    PAL2NAL( ch_cds )
    ch_versions = ch_versions.mix(PAL2NAL.out.versions.first())

    ch_codon_tree = PAL2NAL.out.codon.combine(ch_tree)

    HYPHY_BUSTED( ch_codon_tree )
    HYPHY_RELAX ( ch_codon_tree )
    HYPHY_ABSREL( ch_codon_tree )
    ch_versions = ch_versions.mix(HYPHY_BUSTED.out.versions.first())

    SELECTION_DIAGNOSTICS(
        HYPHY_BUSTED.out.json.map { og, j -> j }.collect(),
        PAL2NAL.out.codon.map { og, c -> c }.collect()
    )

    MK_TEST(
        PAL2NAL.out.codon.map { og, c -> c }.collect(),
        ch_species,
        ch_species
    )

    emit:
    busted      = SELECTION_DIAGNOSTICS.out.summary
    diagnostics = SELECTION_DIAGNOSTICS.out.diagnostics
    mk          = MK_TEST.out.genome_wide
    mk_sweep    = MK_TEST.out.sweep
    versions    = ch_versions
}
