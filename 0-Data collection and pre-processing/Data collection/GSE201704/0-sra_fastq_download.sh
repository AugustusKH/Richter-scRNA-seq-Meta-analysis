#!/bin/bash

# Load sra-tools
module load sra-tools/3.0.3

# Loop through SRR IDs from 18938655 to 18938702
for i in $(seq 18938655 18938702); do
    SRR_ID="SRR${i}"
    echo "Downloading and splitting $SRR_ID..."
    fastq-dump --split-files "$SRR_ID"
done
