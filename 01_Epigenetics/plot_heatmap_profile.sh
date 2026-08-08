#!/bin/bash
## Generic deepTools plotting pipeline
##   For each input group: computeMatrix -> heatmap + profile
##   Group 1: region_sets (reference point: center)
##   Group 2: tss_regions (reference point: TSS, optional)
## Bigwig count and region sets are not limited.
## Usage: fill INPUTS below, then: bash plot_heatmap_profile.sh

set -euo pipefail

## INPUTS ======================================================================
outdir="."
threads=20
flank=3000
refpoint="center"

bws=()              # bigwig files (any number, any marks)
region_sets=()      # region beds for heatmap + profile (any number)
tss_regions=()      # region beds for TSS profile (any number, optional)
## =============================================================================

check_file() {
  local f="$1"
  if [ -z "$f" ]; then echo "Error: empty input path"; return 1; fi
  if [[ "$f" == *\** ]]; then
    compgen -G "$f" >/dev/null && return 0
  elif [ -e "$f" ]; then
    return 0
  fi
  echo "Error: not found: $f"
  return 1
}

validate() {
  local f
  for f in "$@"; do check_file "$f" || exit 1; done
}

plot_both() { # $1: out prefix, $2: reference point, rest: region beds -> matrix + heatmap + profile
  local out="$1" ref="$2"; shift 2
  computeMatrix reference-point \
    -S "${bws[@]}" -R "$@" \
    -a "$flank" -b "$flank" --referencePoint "$ref" \
    -p "$threads" --skipZeros \
    -o "${outdir}/${out}.matrix.gz"
  plotHeatmap -m "${outdir}/${out}.matrix.gz" -out "${outdir}/${out}.pdf" \
    --plotFileFormat pdf --dpi 720 --colorMap YlGnBu --missingDataColor "#FFF6EB"
  plotProfile -m "${outdir}/${out}.matrix.gz" -out "${outdir}/${out}.profile.pdf" --perGroup
}

mkdir -p "$outdir"

[ ${#bws[@]} -gt 0 ] || { echo "Error: bws is empty"; exit 1; }
[ ${#region_sets[@]} -gt 0 ] || { echo "Error: region_sets is empty"; exit 1; }
validate "${bws[@]}" "${region_sets[@]}"

echo "--- Heatmap + profile over region_sets ---"
plot_both regions "$refpoint" "${region_sets[@]}"

if [ ${#tss_regions[@]} -gt 0 ]; then
  echo "--- Heatmap + profile over tss_regions ---"
  validate "${tss_regions[@]}"
  plot_both tss TSS "${tss_regions[@]}"
fi

echo "Done"
