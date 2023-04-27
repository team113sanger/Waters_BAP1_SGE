///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                                 DSL 2                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

nextflow.preview.dsl=2

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                                 USAGE                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

def helpMessage() {
    log.info """
    Usage:

    The typical command for running the pipeline is as follows:

    	nextflow run nf-sge/sge --reads '*_R{1,2}.fastq.gz'


    Mandatory arguments:
      --reads                       Path to input read data (MUST be surrounded with quotes).
      --analysis                    Type of analysis (Options: library or indel)

    Pipeline arguments:
      --no_qc                       Don't run FastQC on raw, trimmed and merged reads.
      --no_trim                     Don't perform adaptor or quality trimming
      --no_merge                    Don't perform paired read merging
      --no_analysis           Don't perform quantification

    Optional arguments:
      --name                        Name for the pipeline run.
      --outdir                      The output directory where the results will be saved.
      --help                        Show pipeline usage.

    For software-specific parameters, please see the relevant configuration files in ./config.
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                            VALIDATE INPUTS                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// Show help message
if ( params.help ) {
	helpMessage()
	exit 0
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = workflow.runName
if (params.name) {
	custom_runName = params.name
}
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)){
	custom_runName = workflow.runName
}

def printErr = System.err.&println

// Validate analysis parameter
def analysis_types = ['library', 'indel']
if ( analysis_types.contains( params.analysis ) == false || !params.analysis ) {
	printErr("Analysis must be provided and can only be one of: " + analysis_types.join(',') + ".")
	exit 1
}

// Validate reads have been provided
if ( !params.reads  ) {
	printErr("Reads must be provided.")
	exit 1
}

// Validate library has been provided
if ( !params.library  ) {
	printErr("Library must be provided.")
	exit 1
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                      Parse input files                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

Channel
    .fromFilePairs( params.reads, checkExists: true )
    .set { raw_reads }

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       	  	  MODULES                         -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

// FastQC
include fastqc as fastqc_raw from './modules/quality_control/common' params(params)
include fastqc as fastqc_trimmed from './modules/quality_control/common' params(params)
include fastqc as fastqc_merged from './modules/quality_control/common' params(params)

// MultiQC
include multiqc from './modules/quality_control/common' params(params)

// read trimming
include './modules/read_trimming/common' params(params)

// read merging
include './modules/read_merging/common' params(params)

// saturation mutagenesis
include './modules/saturation_mutagenesis/common' params(params)

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         SUMMARY INFO                              -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

log.info "${workflow.manifest.name} pipeline - version ${workflow.manifest.version}"
log.info "====================================="

def summary = [:]

summary['Run name']           = custom_runName
summary['Reads']              = params.reads
summary['Results directory']  = params.outdir

summary['Quality control']    = params.no_qc ? 'NO' : 'YES'
summary['Trim reads']         = params.no_trim ? 'NO' : 'YES'
summary['Merge reads']        = params.no_merge ? 'NO' : 'YES'
summary['Library analysis']   = params.analysis == 'library' && !params.no_analysis? 'YES' : 'NO'
summary['Indel analysis']     = params.analysis == 'indel' && !params.no_analysis? 'YES' : 'NO'

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")


///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                           	WORKFLOW                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

workflow qc_trim_and_merge {
	get: raw_reads

	main:
		fastqc_raw( raw_reads )
		(fastqc_raw_zip, fastqc_raw_html) = fastqc_raw.out

		cutadapt( raw_reads )
		(cutadapt_trimmed_reads, cutadapt_logs) = cutadapt.out

		fastqc_trimmed( cutadapt_trimmed_reads )
		(fastqc_trimmed_zip, fastqc_trimmed_html) = fastqc_trimmed.out

		seqprep( cutadapt_trimmed_reads )
		(seqprep_merged_reads) = seqprep.out

		fastqc_merged( seqprep_merged_reads )
		(fastqc_merged_zip, fastqc_merged_html) = fastqc_merged.out

		multiqc_ch = fastqc_raw_zip.mix(cutadapt_logs, fastqc_trimmed_zip, fastqc_merged_zip)
		multiqc( multiqc_ch.collect() )

	emit:
		trimmed_reads = cutadapt_trimmed_reads
		trimmed_read_logs = cutadapt_logs
		merged_reads = seqprep_merged_reads
		multiqc_report = multiqc.out
}

workflow library_analysis {
	get: merged_reads

	main:
		library_quantification(merged_reads)
		(library_counts, library_stats) = library_quantification.out
		library_stats.collectFile(name: 'library_stats_mqc.tsv', newLine: true, keepHeader: true, skip : 1, storeDir: "${workDir}").set{library_stats_merged}

	emit:
		library_counts = library_counts
		library_stats = library_stats_merged
}

workflow indels {
	get: merged_reads

	main:
		indel_alignment(merged_reads, params.library)
		sam_to_cigar_counts(indel_alignment.out)

	emit:
		sam_alignments = indel_alignment.out
		cigar_counts = sam_to_cigar_counts.out
}

workflow {
    main:
        qc_trim_and_merge( raw_reads )
        merged_reads = qc_trim_and_merge.out.merged_reads

        // Library quantification analysis ONLY
        // Processes will be ignored if params.no_analysis or params.analysis is set to indel
        library_analysis( merged_reads )

        // Indel analysis ONLY
        // Processes will be ignored if params.no_analysis or params.analysis is set to library_quantification
        indels( merged_reads )

    publish:
        qc_trim_and_merge.out.trimmed_reads to: "${params.outdir}/cutadapt/trimmed", mode: 'symlink'
        qc_trim_and_merge.out.trimmed_read_logs to: "${params.outdir}/cutadapt/logs", mode: 'copy'
        qc_trim_and_merge.out.merged_reads to: "${params.outdir}/seqprep/merged", mode: 'symlink'
        qc_trim_and_merge.out.multiqc_report to: "${params.outdir}/multiqc", mode: 'copy'
        library_analysis.out.library_counts to: "${params.outdir}/library_quantification", mode: 'copy'
        library_analysis.out.library_stats to: "${params.outdir}/library_quantification", mode: 'copy'
        indels.out.sam_alignments to: "${params.outdir}/indel_analysis", mode: 'symlink'
        indels.out.cigar_counts to: "${params.outdir}/indel_analysis", mode: 'copy'
}
