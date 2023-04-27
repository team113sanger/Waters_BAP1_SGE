# Essentiality phenotyping

## Indel quantification (`nf-sge`)

`nf-sge` is a prototype version of the Nextflow pipeline [QUANTS](https://github.com/cancerit/QUANTS) and was used for indel quantification in the essentiality phenotyping experiment. 

Version 0.0.1 of the `nf-sge` pipeline ([nf-sge__version-0.0.1](nf-sge__version-0.0.1)) was run using [Nextflow](https://www.nextflow.io/) version 19.10.0 and Docker (for dependencies e.g. [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)). Included within this repository is the Dockerfile and scripts used to generate that container [nf-sge__docker-version-0.0.5](nf-sge__docker-version-0.0.5). 

`nf-sge` running[EMBOSS Needle](https://emboss.sourceforge.net/download/) version 6.6.0 was used to globally align reads trimmed with Cutadapt version 2.5 with Python 3.6.8.

To run the pipeline:

```
nextflow [path-to]/nf-sge/main.nf -c experiment.config
```

Below is an example of the `experiment.config` used with `nf-sge` version 0.0.1:

```
docker {
	enabled = true
	runOptions = "-v /home/ubuntu/:/home/ubuntu/"
}

process {
	container = "quay.io/vaofford/nf-sge:0.0.5"
}

params {
	//reads
	reads = [full path to reads]
	library = [full path to library]
	analysis = "indel"

	// cutadapt
	cutadapt_pair_adapters = true
	cutadapt_g = 'file:[full path to R1 adapters]'
	cutadapt_G = 'file:[full path to R2 adapters]'

	cutadapt_error_rate = 0.5
	cutadapt_extra_options = false
}
```

Which requires the adapters in [adapters](adapters) for read trimming with cutadapt.