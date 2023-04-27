#!/usr/bin/env bash

#########################################
#               USAGE                   #
#########################################

usage="
  Convert FASTA to TSV (SeqKit)

  Usage : fa2tsv.sh [options]

  Required:-
    -f <fasta>          : Full path to FASTA file for conversion (expects .fasta or .fa or .txt)
    -i <image path>     : Full path to nf-sge Singularity image [default uses NFSGE_IMAGEDIR and NFSGE_VERSION to build path]
    -h                  : usage

  Optional:-
    -o <output directory>    : Full path to output directory
"

while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts "hf:o:i:" option;
    do
        case "$option" in
                h)  echo "$usage" >&2
                    exit 1
                    ;;
                f)  fasta_file="$OPTARG"
                    ;;
                o)  output_directory="$OPTARG"
                    ;;
                i)  singularity_imagepath="$OPTARG"
                    ;;
                *)  echo "$usage" >&2
                    exit 1
                    ;;

        esac
    done
    shift $((OPTIND-1))
done

if [[ -z "$fasta_file" ]]; then
    echo "FASTA file required (-f)" && exit 1
fi

if [[ ! -f "$fasta_file" ]]; then
    echo "FASTA file does not exist: ${fasta_file}" && exit 1
fi

if [[ -z "$singularity_imagepath" ]]; then
    echo "Singularity imagepath required (-i)" && exit 1
fi

if [[ ! -f "$singularity_imagepath" ]]; then
    echo "Singualrity image does not exist: ${singularity_imagepath}" && exit 1
fi

if [[ -z "$output_directory" ]]; then
    output_directory=$(pwd)
fi

if [[ ! -d "$output_directory" ]]; then
    echo "Output directory does not exist: ${output_directory}" && exit 1
fi

tsv_file=$(echo "${fasta_file}" | sed  -e "s/\.txt$/\.tsv/" -e "s/\.fasta$/\.tsv/" -e "s/\.fa$/\.tsv/")

fa2tsv_cmd="singularity exec ${singularity_imagepath} seqkit fx2tab ${fasta_file} > ${tsv_file}"

eval $fa2tsv_cmd