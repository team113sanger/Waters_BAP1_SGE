///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               SEQPREP                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process seqprep {
    validExitStatus 0
    tag "SeqPrep on $sample_id"

    input:
    tuple sample_id, path(reads)

    output:
    set val(sample_id), path("*.merged.fastq.gz")

    when:
    !params.no_merge

    script:
    prefix = sample_id
    seqprep_cmd = "SeqPrep"

    r1_trimmed_reads = "${sample_id}.R1.seqprep.trimmed.fastq.gz"
    r2_trimmed_reads = "${sample_id}.R2.seqprep.trimmed.fastq.gz"
    r1_discarded_reads = "${sample_id}.R1.discarded.fastq.gz"
    r2_discarded_reads = "${sample_id}.R2.discarded.fastq.gz"
    merged_reads = "${sample_id}.merged.fastq.gz"
    pretty_alignment = "${sample_id}.aln.gz"

    seqprep_cmd = seqprep_cmd + " -f ${reads[0]} -r ${reads[1]}"

    seqprep_cmd = seqprep_cmd + " -1 ${r1_trimmed_reads}"
    seqprep_cmd = seqprep_cmd + " -2 ${r2_trimmed_reads}"
    seqprep_cmd = seqprep_cmd + " -3 ${r1_discarded_reads}"
    seqprep_cmd = seqprep_cmd + " -4 ${r2_discarded_reads}"
    seqprep_cmd = seqprep_cmd + " -s ${merged_reads}"
    seqprep_cmd = seqprep_cmd + " -E ${pretty_alignment}"

    seqprep_cmd = (params.seqprep_extra_options) ? "${seqprep_cmd} ${params.seqprep_extra_options}" : seqprep_cmd    

    """
    $seqprep_cmd
    """
}
