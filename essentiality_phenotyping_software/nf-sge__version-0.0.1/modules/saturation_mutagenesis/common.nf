///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          Library quantification                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process library_quantification {
    validExitStatus 0
    tag "Library quantification: $sample_id"

    input:
    tuple sample_id, path(reads)

    output:
    path "*.merged.counts.tsv"
    path "*.stats.tsv"

    when:
    !params.no_analysis && params.analysis == 'library'

    script:
    quantification_script = "${baseDir}/modules/saturation_mutagenesis/scripts/library_quantification.pl"

    """
    perl "${quantification_script}" -l "${params.library}" -s "${reads[0]}" -w "${params.wildtype_sequence}" -o .
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          Indel alignment                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process indel_alignment {
	validExitStatus 0
	tag "Indel alignment: $sample_id and $library"

	input:
	tuple sample_id, path(reads)
	path(library)

	output:
	set val(sample_id), path("*.sam")

	when:
	!params.no_analysis && params.analysis == 'indel'

	script:
	needle_cmd = "needle -auto"
	needle_cmd = needle_cmd + " -asequence ${library}"
	needle_cmd = needle_cmd + " -bsequence reads.fq"
	needle_cmd = needle_cmd + " -sformat2 ${params.needle_sformat2}"
	needle_cmd = needle_cmd + " -aformat3 ${params.needle_aformat3}"
	needle_cmd = needle_cmd + " -aname3 ${sample_id}"
	needle_cmd = needle_cmd + " -aextension3 ${params.needle_aextension3}"
	needle_cmd = needle_cmd + " -gapopen ${params.needle_gapopen}"
	needle_cmd = needle_cmd + " -gapextend ${params.needle_gapextend}"
	needle_cmd = needle_cmd + " -endopen ${params.needle_endopen}"
	needle_cmd = needle_cmd + " -endextend ${params.needle_endextend}"
	needle_cmd = (params.needle_endweight) ? "${needle_cmd} -endweight" : needle_cmd

	"""
	zcat ${reads} > reads.fq
	$needle_cmd
	rm reads.fq
	"""
}

process sam_to_cigar_counts {
	validExitStatus 0
	tag "SAM to CIGAR counts: $sample_id"

	input:
	tuple sample_id, path(sam)

	output:
	set val(sample_id), path("*.counts.tsv")

	when:
	!params.no_analysis && params.analysis == 'indel'

	script:
	"""
	awk -F"\t" '{print \$6}' ${sam} | sort | uniq -c | sort -nr > ${sample_id}.counts.tsv
	"""
}
