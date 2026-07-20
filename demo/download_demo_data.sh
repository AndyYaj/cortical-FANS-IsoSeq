#!/bin/bash

set -e

mkdir -p data

wget -O data/276_GABA_demo.subreads.bam \
"https://zenodo.org/records/20301473/files/276_GABA_demo.subreads.bam?download=1"

wget -O data/276_GABA_demo.subreads.bam.pbi \
"https://zenodo.org/records/20301473/files/276_GABA_demo.subreads.bam.pbi?download=1"
