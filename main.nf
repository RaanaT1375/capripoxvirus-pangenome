#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    capripoxvirus-pangenome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Pangenome, recombination and selection analysis of the Capripoxvirus genus.

    Github : https://github.com/RaanaT1375/capripoxvirus-pangenome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { READ_PREPROCESSING } from './subworkflows/local/read_preprocessing'
include { PANGENOME          } from './subworkflows/local/pangenome'
include { PHYLOGENOMICS      } from './subworkflows/local/phylogenomics'
include { POPULATION         } from './subworkflows/local/population'
include { SELECTION          } from './subworkflows/local/selection'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP AND PARAMETER VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
    log.info """
    ==========================================================================
     capripoxvirus-pangenome  v${workflow.manifest.version}
    ==========================================================================

    Usage:

      nextflow run RaanaT1375/capripoxvirus-pangenome \\
          --input samplesheet.csv \\
          --outdir results \\
          -profile slurm,conda

    Required:
      --input          Samplesheet CSV. Columns:
                         sample,species,host,country,year,condition,
                         assembly,fastq_1,fastq_2
                       Provide either 'assembly' (a FASTA path) or the FASTQ
                       columns for each row, not both.
      --outdir         Output directory.

    References:
      --lsdv_reference     LSDV reference genome FASTA   (default AF325528.1)
      --lsdv_proteome      LSDV reference proteome FAA   (GCF_000839805.1)
      --host_bowtie2_index Directory of prebuilt Bowtie2 host indices
      --bacteria_mmi       minimap2 index of bacterial reference genomes

    Filtering:
      --min_contig_len         500     bp, first-pass assemblies
      --min_contig_len_rescue  1000    bp, bacterially depleted re-assemblies
      --cds_flag_threshold     200     predicted CDS above which an assembly is
                                       flagged for bacterial depletion
      --contamination_max      0.20    max fraction of unmapped proteins
      --max_missing_near_core  9       max missing near-core genes

    Analysis:
      --core_threshold      0.99   presence fraction defining the strict core
      --gubbins_iterations  5
      --scoary_permutations 1000
      --mk_freq_cutoff      0.05   minor-allele frequency cutoff
      --mk_freq_sweep       '0,0.01,0.02,0.05,0.10,0.15'
      --skip_selection      false
      --skip_phylogeography false

    Profiles: conda, mamba, singularity, docker, slurm, test
    ==========================================================================
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

if (!params.input)  { exit 1, "ERROR: --input samplesheet not specified. Use --help." }
if (!params.outdir) { exit 1, "ERROR: --outdir not specified. Use --help." }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMPLESHEET PARSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parseSamplesheet(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample) exit 1, "ERROR: samplesheet row without a 'sample' value."

            def meta = [
                id       : row.sample,
                species  : row.species   ?: 'NA',
                host     : row.host      ?: 'NA',
                country  : row.country   ?: 'NA',
                year     : row.year      ?: 'NA',
                condition: row.condition ?: 'field'
            ]

            def hasAssembly = row.assembly?.trim()
            def hasReads    = row.fastq_1?.trim()

            if (hasAssembly && hasReads)
                exit 1, "ERROR: ${row.sample} declares both an assembly and FASTQ files. Choose one."
            if (!hasAssembly && !hasReads)
                exit 1, "ERROR: ${row.sample} declares neither an assembly nor FASTQ files."

            if (hasAssembly) {
                return [ 'assembly', meta, file(row.assembly, checkIfExists: true) ]
            }
            meta.single_end = !row.fastq_2?.trim()
            def reads = meta.single_end
                ? [ file(row.fastq_1, checkIfExists: true) ]
                : [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
            return [ 'reads', meta, reads ]
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAPRIPOX {

    ch_versions = Channel.empty()

    // ---- inputs -----------------------------------------------------------
    ch_rows = parseSamplesheet(params.input)

    ch_reads = ch_rows.filter { it[0] == 'reads'    }.map { [ it[1], it[2] ] }
    ch_given = ch_rows.filter { it[0] == 'assembly' }.map { [ it[1], it[2] ] }

    ch_lsdv_ref      = Channel.fromPath(params.lsdv_reference, checkIfExists: true).collect()
    ch_lsdv_proteome = Channel.fromPath(params.lsdv_proteome,  checkIfExists: true).collect()

    // ---- 01-04  reads to annotated assemblies -----------------------------
    READ_PREPROCESSING( ch_reads, ch_lsdv_proteome )
    ch_versions = ch_versions.mix(READ_PREPROCESSING.out.versions)

    ch_assemblies = ch_given.mix(READ_PREPROCESSING.out.assemblies)

    // ---- 05-06  QC and pangenome ------------------------------------------
    PANGENOME( ch_assemblies, ch_lsdv_proteome )
    ch_versions = ch_versions.mix(PANGENOME.out.versions)

    // ---- 07-08  supermatrix and ML tree -----------------------------------
    PHYLOGENOMICS( PANGENOME.out.single_copy_orthogroups, PANGENOME.out.faa_clean )
    ch_versions = ch_versions.mix(PHYLOGENOMICS.out.versions)

    // ---- 09-13  recombination, structure, phylogeography, pan-GWAS --------
    POPULATION(
        PANGENOME.out.fna_clean,
        PHYLOGENOMICS.out.tree,
        PANGENOME.out.presence_absence,
        ch_lsdv_ref
    )
    ch_versions = ch_versions.mix(POPULATION.out.versions)

    // ---- 14-16  selection and population genetics -------------------------
    if (!params.skip_selection) {
        SELECTION(
            PHYLOGENOMICS.out.gene_alignments,
            PHYLOGENOMICS.out.tree,
            PANGENOME.out.species_assignment
        )
        ch_versions = ch_versions.mix(SELECTION.out.versions)
    }

    // ---- provenance -------------------------------------------------------
    ch_versions
        .unique()
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")
}

workflow {
    CAPRIPOX()
}

workflow.onComplete {
    log.info """
    --------------------------------------------------------------------------
     Completed  : ${workflow.success ? 'OK' : 'FAILED'}
     Duration   : ${workflow.duration}
     Results    : ${params.outdir}
     Report     : ${params.outdir}/pipeline_info/execution_report.html
    --------------------------------------------------------------------------
    """.stripIndent()
}
