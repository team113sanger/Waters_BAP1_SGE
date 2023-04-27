#!/usr/bin/env bash

#########################################
#               USAGE                   #
#########################################

usage="
  Rename paired MiSeq FASTQ files using tab-delimited mapping file

  Usage :  renameMiSeqFastq.sh [options] -i <mapping file>

  Required:-
    -i <mapping file>      : tab-delimited mapping file 
    -h                     : usage

  Optional:-
    -f <FASTQ directory>   : directory where original FASTQ files are stored [Default: .]
    -o <output directory>  : directory where renamed FASTQ files will be saved [Default: .]
    -e <file extension>    : file extension [Default: .fastq.gz]
"

while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts "hi:f:o:e:" option;
    do
        case "$option" in
                h)  echo "$usage" >&2
                    exit 1
                    ;;
				f)	fq_directory="$OPTARG"
					;;
				o)  output_directory="$OPTARG"
                    ;;
				e)	file_ext="$OPTARG"
					;;
                i)  id_file="$OPTARG"
                    ;;
                *)  echo "$usage" >&2
                    exit 1
                    ;;

        esac
    done
    shift $((OPTIND-1))
done

if [[ -z "$id_file" ]]; then
    echo "Tab-delimited mapping file required (-i)" && exit 1
fi

if [[ ! -f "$id_file" ]]; then
    echo "Tab-delimited mapping file does not exist: ${id_file}" && exit 1
fi

ncol=$(awk -F'\t' 'NR==1{print NF}' $id_file )
if [[ ncol -lt 2 ]]; then
	echo "Tab-delimited mapping file has less than 2 columns: ${id_file}" && exit 1
elif [[ ncol -gt 2 ]]; then
	echo "Tab-delimited mapping file has more than 3 columns, only taking the first two: ${id_file}"
fi

if [[ -z "$fq_directory" ]]; then
	echo "$PWD"
	fq_directory=$PWD
    echo "No FASTQ directory given, searching ${fq_directory}"
fi

if [[ ! -d "$fq_directory" ]]; then
	echo "FASTQ directory does not exist: ${fq_directory}" && exit 1
fi

if [[ -z "$output_directory" ]]; then
	output_directory=$PWD
	echo "No output directory given, renaming files to ${output_directory}"
fi

if [[ ! -d "$output_directory" ]]; then
    echo "Output directory does not exist: ${output_directory}" && exit 1
fi

if [[ -z "$file_ext" ]]; then
	file_ext=".fastq.gz"
fi

rename_log="${output_directory}/rename_fq.log"

while read line; do
	old_prefix=$(echo "$line" | cut -d $'\t' -f 1)
	new_prefix=$(echo "$line" | cut -d $'\t' -f 2)

	if [[ -z $old_prefix ]] || [[ -z $new_prefix ]]; then
		echo "Could not extract old and new prefices from tab-delimited mapping file: ${id_file}" && exit 1
	fi

	targets=($( find ${fq_directory} -name "${old_prefix}*${file_ext}" ))
	
	if [[ ${#targets[@]} -lt 1 ]]; then
		echo "Could not find any files matching ${old_prefix}*${file_ext} in ${fq_directory}" >> $rename_log
	fi

	for old_path in "${targets[@]}"
	do
		old_fname=$(basename $old_path)
		new_fname="${old_fname/$old_prefix/$new_prefix}"
		new_fname="${new_fname/$file_ext/.fastq.gz}"
		new_path="${output_directory}/${new_fname}"

		echo "Moving ${old_path} to ${new_path}" >> $rename_log
		mv_cmd="mv ${old_path} ${new_path}"		
		mv_output="$($mv_cmd)"
	done
	
done < $id_file

