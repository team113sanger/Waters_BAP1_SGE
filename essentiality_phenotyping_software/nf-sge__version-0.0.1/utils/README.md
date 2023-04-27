# sge_pipeline


## Downloading project data from BaseSpace

Before you can download your data using the BaseSpace Sequence Hub CLI using the provided scripts, you will need to have done the following:

1) Download BaseSpace Sequence Hub CLI (**bs**) by following the instructions [here](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview).
2) Make sure that **bs** is in your PATH by running `bs --version`
3) Authenticate your account by running `bs authenticate`. This will prompt you to go to a validation URL and then return your account details.
4) Check that the authentication worked by running `bs whoami`. 
5) Find the project name whose data you want to download from those given by `bs list projects`.

To download FASTQ files for a BaseSpace project you can use `scripts/download_basespace_fastq.sh`. 

```
scripts/download_basespace_fastq.sh -h
```

Example command:

```
scripts/download_basespace_fastq.sh -p <project_name> -o <output_directory>
```

The script will create a project folder in the output directory (e.g. project_name). It will then download a compressed file (e.g. project_name.tar.gz) in that directory which contains the gzipped FASTQ files for the project into the project folder.  Finally, the script will extract the gzipped FASTQ files into the project directory (e.g. output_directory/project_name/reads1.fastq.gz). 

SGE analysis scripts will assume that all of the FASTQ files you want to analyse are in a single folder.


## Getting the software

If you have Docker installed, the software used in this analysis is available as a Docker container:

```
docker pull quay.io/vaofford/nf-sge:latest
```

If you have Singularity installed, you can also generate a Singularity image from the Docker container with:

```
singularity pull docker://quay.io/vaofford/nf-sge:latest
```

Docker is not available on the Sanger compute farm. However, it can be used in the Sanger FCE. For users of the Sanger farm, we have included a script to generate the Singularity image: 

```
scripts/docker2sge.sh 
```

If you want a particular tagged version, you can set `NFSGE_VERSION`. To set the output directory for the Singularity image, you can set `NFSGE_IMAGEDIR`. 

Alternatively, it is possible to install all of the software components from source. That can be quite painful, but, please, fire ahead if that's your preferred method!  


## Renaming FASTQ files

Sometimes files come with prefices or suffices which aren't very informational (e.g. 1, 2, 3.....). To rename FASTQ files, you will need a tab-delimited file with the expected sample prefix in the first column and the prefix that you want to replace it with in the second column. For example, to rename *1_S1_L001_R2_001.fastq.gz* to *day_5_H03_BRCA1_BRCA1_ex_9_R2.fastq.gz* you would need to have *1_S1_L001* in the first column and *day_5_H03_BRCA1_BRCA1_ex_9* in the second column.

To rename FASTQ files:

```
scripts/renameMiSeqFastq.sh -i <mapping file>
```

The script will search and rename files to your current working directory unless you specify a FASTQ directory (`-f`) and/or an output directory (`-o`) respectively. Files suffices are replaced with ".fastq.gz" using the file extension (`-e`) which is set to "001.fastq.gz".  

## Converting a library from FASTA to TSV

The pipeline assumes that the library is a tab-delimited file where the first column is the sequence header and the second column is the sequence. A FASTA formatted library can be converted using `fa2tsv.sh` which will run SeqKit.  This will replace the original file extension with ".tsv". For example, if the FASTA input file was called _library.fa_ the TSV output file will be called _library.tsv_.

To convert a FASTA library to TSV:

```
scripts/fa2tsv.sh -f <library.fa> -i <singularity_image_path>
```

_N.B. your library file must have one of the following file extensions: .fa, .fasta or .txt._
