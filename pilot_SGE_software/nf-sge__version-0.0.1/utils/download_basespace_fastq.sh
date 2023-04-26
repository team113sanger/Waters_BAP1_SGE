#!/usr/bin/env bash

#########################################
#               USAGE                   #
#########################################

usage="
Download project FASTQ files from BaseSpace

  Usage : download_basespace_fastq.sh [options] 

  Required:-
    -p <BaseSpace project name>      : BaseSpace project name
    -h                               : usage
    
  Optional:-
    -s <BaseSpace API hostname>      : BaseSpace API hostname (server)
    -a <BaseSpace API access token>  : Manual BaseSpace API access token
    -o <output directory>            : Downloaded FASTQ destination directory [default .]
"

while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts "hp:s:a:o:" option;
    do
        case "$option" in
                h)  echo "$usage" >&2
                    exit 1
                    ;;
                p)  bs_project_name="$OPTARG"
                    ;;
                s)  bs_hostname="$OPTARG"
                    ;;
                a)  bs_access_token="$OPTARG"
                    ;;
                o)  output_directory="$OPTARG"
                    ;;
                *)  echo "$usage" >&2
                    exit 1
                    ;;

        esac
    done
    shift $((OPTIND-1))
done

if [[ ! $(type -P "bs") ]]; then 
    echo "bs not in PATH" && exit 1
fi

if [[ -z $(bs whoami) ]]; then
    echo "bs not authenicated" && exit 1
fi

if [[ -z $(bs list projects --) ]]; then
    echo "bs not authenicated" && exit 1
fi

if [[ -z "$output_directory" ]]; then
    output_directory=$(pwd)
fi

if [[ ! -d "$output_directory" ]]; then
    echo "Output directory does not exist: ${output_directory}" && exit 1
fi

if [[ -z "$bs_project_name" ]]; then 
    echo "Project name missing" && exit 1
fi

project_fastq_directory="${output_directory}/${bs_project_name}"

if [[ -d "$project_fastq_directory" ]]; then
    echo "Project FASTQ directory already exists: ${project_fastq_directory}" && exit 1
else
    echo "Creating project directory: ${project_fastq_directory}"
    mkdir $project_fastq_directory
fi

if [[ ! -d "$project_fastq_directory" ]]; then
    echo "Project FASTQ does not exist: ${project_fastq_directory}" && exit 1
fi

original_directory=$(pwd)
cd $project_fastq_directory

bs_extenstion=".fastq.gz"
bs_dnld_cmd="bs download project --name=${bs_project_name} --extension=${bs_extenstion} -z --output=${bs_project_name}"

if [[ ! -z "$bs_hostname" ]]; then
    bs_dnld_cmd="${bs_dnld_cmd} --api-server=${bs_hostname}"
fi

if [[ ! -z "$bs_access_token" ]]; then
    bs_dnld_cmd="${bs_dnld_cmd} --access-token=${bs_access_token}"
fi

echo "Running BaseSpace CLI download command: ${bs_dnld_cmd}"
dnld_log="${output_directory}/${bs_project_name}_bs_dnld.log"
bs_dnld_output="$($bs_dnld_cmd > $dnld_log)"

if [[ ! -z "$bs_dnld_output" ]]; then
    echo "DOWNLOAD FAILED (command error): ${bs_dnld_output}" && exit 1
fi

bs_dnld_fname="${bs_project_name}.tar.gz"
if [[ ! -e "${bs_dnld_fname}" ]]; then
    echo "DOWNLOAD FAILED (file does not exist): ${bs_dnld_fname}" && exit 1
fi

echo "Uncompressing downloaded data: ${bs_dnld_fname} to ${project_fastq_directory}"
tar_cmd="tar -xf ${bs_dnld_fname} --strip-components=2"
tar_log="${output_directory}/${bs_project_name}_tar.log"
tar_output="$($tar_cmd > $tar_log)"

if [[ ! -z "$bs_dnld_fname" ]]; then
    echo "Removing tarball: ${bs_dnld_fname}"
    rm ${bs_dnld_fname}
fi

cd $original_directory

echo "SUCCESS: FASTQ files (.fastq.gz) have been downloaded to ${project_fastq_directory}"
