///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               FASTQC                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process fastqc {
    validExitStatus 0,143
    tag "FastQC on $sample_id"

    input:
    tuple sample_id, path(reads)

    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    when:
    !params.no_qc

    script:
    """
    fastqc --quiet $reads
    """
}

///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                               MULTIQC                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////

process multiqc {
    validExitStatus 0
    tag "MultiQC"

    input:
    path("*")

    output:
    path ("*.html")

    when:
    !params.no_qc

    multiqc_html = "${workflow.manifest.name}.multiqc.html"
    multiqc_html = multiqc_html.replaceAll(/ /, "_")

    script:

    if( !params.multiqc_config ) 
        """
        multiqc -i "${workflow.manifest.name}" -n "${multiqc_html}" . 
        """
    else
        """
        multiqc -i "${workflow.manifest.name}" -n "${multiqc_html}" -c "${params.multiqc_config}" . 
        """
}
