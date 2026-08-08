#!/bin/bash
## HOMER motif finding + peak annotation
##   findMotifsGenome.pl: de novo motif discovery in peaks
##   annotatePeaks.pl:    peak-to-gene annotation
## Usage: $0 <sample> [peaks.bed]
##   peaks.bed: BED file, default homer_peaks.tmp

set -euo pipefail

sample="${1:?Usage: $0 <sample> [peaks.bed]}"
peaks="${2:-homer_peaks.tmp}"

genome="hg38"
motif_len="8,10,12"

[ -f "$peaks" ] || { echo "Error: peaks file not found: $peaks"; exit 1; }
for cmd in findMotifsGenome.pl annotatePeaks.pl; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found (HOMER not in PATH)"; exit 1; }
done

findMotifsGenome.pl "$peaks" "$genome" "${sample}_motifDir" -len "$motif_len"

annotatePeaks.pl "$peaks" "$genome" 1> "${sample}.peakAnn.xls" 2> "${sample}.annLog.txt"

echo "Done"
