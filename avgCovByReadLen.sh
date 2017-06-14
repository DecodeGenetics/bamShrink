#!/usr/bin/bash
set -e
set -o pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: bash avgCovByReadLen.sh <BAM>"
  exit 1
fi

bam=$1 # Path to the target BAM file

if [[ ! -f $bam ]]; then
  echo "ERROR: Could not find BAM at path '${bam}'"
  exit 1
fi

if [[ ! -f ${bam}.bai ]]; then
  echo "ERROR: BAM file needs to be indexed first."
  exit 1
fi

samtools idxstats $bam | head -n -1 | awk '{sum+=$3+$4; ref+=$2;} END{print sum/ref}'
