#!/bin/sh

usage() {  # Print a help message.
  echo "Usage: $0 [ -p pipeline_version ] [ -n nextflow_version ] [ -b basespace_version ] [-i install_dir] [-o openstack]" 1>&2
}

exit_error() {
  usage
  exit 1
}

pipeline_version="v1.0.0"
nextflow_version="19.10.0"
basespace_version="1.1.0"
install_dir="/home/ubuntu/bin"
openstack=1

while getopts ":p:n:b:o:i:h" option; do
case "${option}"
in
	p) pipeline_version=${OPTARG};;
	n) nextflow_version=${OPTARG};;
  b) basespace_version=${OPTARG};;
  o) openstack=${OPTARG};;
  i) install_dir=${OPTARG};;
  h) usage;;
esac
done

# Install Java version 8 if setting up on openstack
if [ $openstack -eq 1 ]
then
  apt-get update
  apt-get install -y openjdk-8-jre-headless
fi

# Make bin directory for executables.
mkdir ${install_dir}

# Add bin directory to PATH so executables are available.
export PATH=$PATH:$install_dir

# Download and install nextflow.
wget "http://www.nextflow.io/releases/v${nextflow_version}/nextflow-${nextflow_version}-all" -O ${install_dir}/nextflow

# Download and install Basespace CLI.
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/${basespace_version}/amd64-linux/bs?bt_package=latest" -O ${install_dir}/bs

# Update bin directory permissions.
chmod 777 ${install_dir}/*

# Clone SGE pipeline repository (need to have SSH keys set up or use HTTPS command instead). Version 1.0.0 was used for this analysis (commit 4bb3116d).
git clone --branch ${pipeline_version} --single-branch https://vo1:D7bp1LfdBmWj_cVfeFd4@gitlab.internal.sanger.ac.uk/nextflow_pipelines/nf-sge.git ${install_dir}/nf-sge
