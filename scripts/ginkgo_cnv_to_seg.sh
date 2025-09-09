#!/bin/bash
# Usage
#
#   bash path/to/ginkgo_cnv_to_seg.sh P33410_10*CNV1.tsv > Gingko_CNV1.seg

set -euo pipefail

# Header for IGV
echo "#type=COPY_NUMBER"

for file in "$@"
do
    n=2
    # If first file, print header
    if [ "$file" == "$1" ]; then
        n=1
    fi

    tail -n+$n $file | awk '{OFS="\t"; print $4, $1, $2, $3, $5}'

done
