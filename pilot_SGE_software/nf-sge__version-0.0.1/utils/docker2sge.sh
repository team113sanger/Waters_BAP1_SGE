#!/usr/bin/env bash

container_version="latest"
image_version=$(date '+%d%m%Y%H%M%S');

if [[ ! -v NFSGE_VERSION ]] || [[ -z "$NFSGE_VERSION" ]]; then
    echo "NFSGE_VERSION is not set or empty, getting latest Docker container"
else
    container_version=$NFSGE_VERSION
    image_version=$NFSGE_VERSION
fi

if [[ ! -v NFSGE_IMAGEDIR ]] || [[ -z "$NFSGE_IMAGEDIR" ]]; then
    image_dest="${PWD}/nf-sge_${image_version}.sif"
else
    image_dest="${NFSGE_IMAGEDIR}/nf-sge_${image_version}.sif"
fi

container_uri="docker://quay.io/vaofford/nf-sge:${container_version}"

echo "Downloading Docker container ${container_uri} to ${image_dest}"
module load ISG/singularity/03.2.0
singularity pull $image_dest $container_uri
module unload ISG/singularity

if [ -f "$image_dest" ]; then
    echo "Singularity container successfully built: ${image_dest}"
else
    echo "ERROR: Singularity container build failed"
fi