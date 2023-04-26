///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                              CUTADAPT                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process cutadapt {
  validExitStatus 0
  tag "cutadapt on $sample_id"

  input:
    tuple sample_id, path(reads) 

  output:
    set val(sample_id), path("*.trimmed.fastq.gz")
    path ("*.cutadapt.out")

  when:
    !params.no_trim

	script:
    prefix = sample_id
    cutadapt_cmd = "cutadapt"

    cutadapt_cmd = (params.cutadapt_cores) ? "${cutadapt_cmd} --cores ${params.cutadapt_cores}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_error_rate) ? "${cutadapt_cmd} --error-rate ${params.cutadapt_error_rate}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_times) ? "${cutadapt_cmd} --times ${params.cutadapt_times}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_overlap) ? "${cutadapt_cmd} --overlap ${params.cutadapt_overlap}" : cutadapt_cmd

    cutadapt_cmd = (params.cutadapt_a) ? "${cutadapt_cmd} -a ${params.cutadapt_a}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_g) ? "${cutadapt_cmd} -g ${params.cutadapt_g}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_b) ? "${cutadapt_cmd} -b ${params.cutadapt_b}" : cutadapt_cmd

    cutadapt_cmd = (params.cutadapt_pair_adapters) ? "${cutadapt_cmd} --pair-adapters" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_pair_filter) ? "${cutadapt_cmd} --pair-filter ${params.cutadapt_pair_filter}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_A) ? "${cutadapt_cmd} -A ${params.cutadapt_A}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_G) ? "${cutadapt_cmd} -G ${params.cutadapt_G}" : cutadapt_cmd
    cutadapt_cmd = (params.cutadapt_B) ? "${cutadapt_cmd} -B ${params.cutadapt_B}" : cutadapt_cmd

    cutadapt_cmd = (params.cutadapt_extra_options) ? "${cutadapt_cmd} ${params.cutadapt_extra_options}" : cutadapt_cmd

    cutadapt_trimmed_r1 = sample_id + '.R1.trimmed.fastq.gz'
    cutadapt_trimmed_r2 = sample_id + '.R2.trimmed.fastq.gz'

    cutadapt_cmd = cutadapt_cmd + " --output ${cutadapt_trimmed_r1}"
    cutadapt_cmd = cutadapt_cmd + " --paired-output ${cutadapt_trimmed_r2}"

    cutadapt_log = prefix + ".cutadapt.out"

    """
    $cutadapt_cmd ${reads[0]} ${reads[1]} > $cutadapt_log
    """
}