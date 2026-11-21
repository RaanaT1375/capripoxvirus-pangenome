//
// Stages 01-04: raw reads to annotated assemblies.
//
// Two-pass design. Every dataset is trimmed, host-depleted and assembled with
// SPAdes --careful, then annotated. Assemblies whose predicted CDS count
// exceeds the expected capripoxvirus repertoire are flagged, their reads are
// depleted of bacterial sequence with minimap2, and they are re-assembled with
// --only-assembler and re-annotated. A blanket read-level decontamination is
// deliberately avoided because it would discard lineage-specific accessory
// genes, which are the compartment of interest.
//

include { FASTP                          } from '../../modules/local/fastp'
include { BOWTIE2_HOST_DEPLETION         } from '../../modules/local/bowtie2_host_depletion'
include { MINIMAP2_BACTERIAL_DEPLETION   } from '../../modules/local/minimap2_bacterial_depletion'
include { SPADES as SPADES_FIRST         } from '../../modules/local/spades'
include { SPADES as SPADES_RESCUE        } from '../../modules/local/spades'
include { PROKKA as PROKKA_FIRST         } from '../../modules/local/prokka'
include { PROKKA as PROKKA_RESCUE        } from '../../modules/local/prokka'
include { CDS_SCREEN                     } from '../../modules/local/cds_screen'

workflow READ_PREPROCESSING {

    take:
    ch_reads          // [ meta, [ fastq ] ]
    ch_proteome       // reference proteome for Prokka

    main:
    ch_versions = Channel.empty()

    FASTP( ch_reads )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    ch_host_index = params.host_bowtie2_index
        ? Channel.fromPath(params.host_bowtie2_index, checkIfExists: true).collect()
        : Channel.value([])

    BOWTIE2_HOST_DEPLETION( FASTP.out.reads, ch_host_index )
    ch_versions = ch_versions.mix(BOWTIE2_HOST_DEPLETION.out.versions.first())

    // ---- first pass --------------------------------------------------------
    SPADES_FIRST( BOWTIE2_HOST_DEPLETION.out.reads, 'careful' )
    PROKKA_FIRST( SPADES_FIRST.out.contigs, ch_proteome )
    CDS_SCREEN( PROKKA_FIRST.out.tsv )

    ch_flagged = CDS_SCREEN.out.verdict
        .filter { meta, flag, tsv -> flag == 'flagged' }
        .map    { meta, flag, tsv -> meta.id }

    ch_clean_first = PROKKA_FIRST.out.faa
        .join( CDS_SCREEN.out.verdict.map { meta, flag, tsv -> [ meta, flag ] } )
        .filter { meta, faa, flag -> flag == 'clean' }
        .map    { meta, faa, flag -> [ meta, faa ] }

    // ---- second pass, flagged datasets only --------------------------------
    if (params.bacteria_mmi) {
        ch_bact = Channel.fromPath(params.bacteria_mmi, checkIfExists: true).collect()

        ch_to_rescue = BOWTIE2_HOST_DEPLETION.out.reads
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .join( ch_flagged.map { id -> [ id, true ] } )
            .map { id, meta, reads, flagged -> [ meta, reads ] }

        MINIMAP2_BACTERIAL_DEPLETION( ch_to_rescue, ch_bact )
        SPADES_RESCUE( MINIMAP2_BACTERIAL_DEPLETION.out.reads, 'rescue' )
        PROKKA_RESCUE( SPADES_RESCUE.out.contigs, ch_proteome )

        ch_versions = ch_versions.mix(MINIMAP2_BACTERIAL_DEPLETION.out.versions.first())
        ch_faa = ch_clean_first.mix( PROKKA_RESCUE.out.faa )
        ch_ffn = PROKKA_FIRST.out.ffn.mix( PROKKA_RESCUE.out.ffn )
    } else {
        log.warn "No --bacteria_mmi provided: flagged assemblies are reported but not re-processed."
        ch_faa = PROKKA_FIRST.out.faa
        ch_ffn = PROKKA_FIRST.out.ffn
    }

    ch_versions = ch_versions.mix(SPADES_FIRST.out.versions.first())
    ch_versions = ch_versions.mix(PROKKA_FIRST.out.versions.first())

    emit:
    assemblies = ch_faa.join(ch_ffn, remainder: true).map { meta, faa, ffn -> [ meta, faa ] }
    faa        = ch_faa
    ffn        = ch_ffn
    cds_counts = CDS_SCREEN.out.counts
    versions   = ch_versions
}
