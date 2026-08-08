#!/bin/bash
## ATAC-seq paired-end pipeline: mapping -> dedup -> bigwig -> peak calling
## Usage: $0 [sample.list]
##   sample.list: <sample> <fastq1> <fastq2> per line

set -euo pipefail

## INPUTS ======================================================================
sample_list="${1:-./sample.list}"

bowtie2_index=""          # bowtie2 index prefix (full path)
atac_included=""          # restriction regions BED

threads_bowtie2=20
threads_samtools=15
threads_light=4
threads_index=10
threads_bamcoverage=20

max_fragment=2000
multimapping=10
min_mapq=30
filter_flag=524
dup_flag=1024
java_mem="-Xmx64G -Xms8G"

env_gatk="gatk"
env_chipseq="chipseq"

genome="hs"
peak_q=0.01
macs_shift=-100
macs_extsize=200

outdir="mapping"
peak_dir="peak_calling"
## =============================================================================

bw_dir="${outdir}/bw"

[ -f "$sample_list" ] || { echo "Error: sample list not found: $sample_list"; exit 1; }
[ -n "$bowtie2_index" ] || { echo "Error: bowtie2_index is empty"; exit 1; }
[ -f "$atac_included" ] || { echo "Error: atac_included not found: $atac_included"; exit 1; }
ls "${bowtie2_index}".*.bt2 >/dev/null 2>&1 || { echo "Error: bowtie2 index not found: $bowtie2_index"; exit 1; }

for cmd in bowtie2 samtools conda; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH"; exit 1; }
done

activate() { # activate conda env if specified
  local env="$1"
  [ -n "$env" ] || return 0
  if ! type __conda_activate >/dev/null 2>&1; then
    local conda_bin
    conda_bin="$(command -v conda)"
    source "$(dirname "$(dirname "$conda_bin")")/etc/profile.d/conda.sh"
  fi
  conda activate "$env"
}

map_track_a() { # multi-mapping (-k) alignment, name/coordinate sorted, dedup
  local sample="$1" fq1="$2" fq2="$3"
  local sam="${outdir}/${sample}.sam"
  bowtie2 --very-sensitive -k "$multimapping" -X "$max_fragment" \
    --threads "$threads_bowtie2" -x "$bowtie2_index" -1 "$fq1" -2 "$fq2" -S "$sam"
  samtools view --threads "$threads_samtools" -F "$filter_flag" -L "$atac_included" -b "$sam" \
    | samtools sort -n --threads "$threads_light" -O BAM -o "${outdir}/${sample}_nSorted.bam"
  samtools sort --threads "$threads_samtools" -O BAM \
    -o "${outdir}/${sample}_pSorted.bam" "${outdir}/${sample}_nSorted.bam"
  samtools index -@ "$threads_index" "${outdir}/${sample}_pSorted.bam"
  activate "$env_gatk"
  gatk --java-options "$java_mem" MarkDuplicates \
    --INPUT "${outdir}/${sample}_pSorted.bam" \
    --OUTPUT "${outdir}/${sample}_MD.bam" \
    --METRICS_FILE "${outdir}/${sample}_MD.txt"
  samtools index -@ "$threads_index" "${outdir}/${sample}_MD.bam"
  rm -f "$sam"
}

map_track_b() { # unique alignment + dedup removal
  local sample="$1" fq1="$2" fq2="$3"
  local sam="${outdir}/${sample}.nk.sam"
  bowtie2 --very-sensitive -X "$max_fragment" \
    --threads "$threads_bowtie2" -x "$bowtie2_index" -1 "$fq1" -2 "$fq2" -S "$sam"
  samtools view --threads "$threads_light" -F "$filter_flag" -q "$min_mapq" -L "$atac_included" -b "$sam" \
    | samtools sort --threads "$threads_light" -O BAM -o "${outdir}/${sample}_Sorted.nk.bam"
  activate "$env_gatk"
  gatk --java-options "$java_mem" MarkDuplicates \
    --INPUT "${outdir}/${sample}_Sorted.nk.bam" \
    --OUTPUT "${outdir}/${sample}_MD.nk.bam" \
    --METRICS_FILE "${outdir}/${sample}_md.metrics"
  samtools index -@ "$threads_index" "${outdir}/${sample}_MD.nk.bam"
  samtools view --threads "$threads_samtools" -F "$dup_flag" -b \
    -o "${outdir}/${sample}_RD.nk.bam" "${outdir}/${sample}_MD.nk.bam"
  samtools index -@ "$threads_index" "${outdir}/${sample}_RD.nk.bam"
  rm -f "$sam"
}

to_bigwig() {
  local sample="$1"
  activate "$env_chipseq"
  bamCoverage --numberOfProcessors "$threads_bamcoverage" --normalizeUsing RPKM \
    --outFileFormat bigwig -b "${outdir}/${sample}_RD.nk.bam" -o "${bw_dir}/${sample}_RPKM.bw"
}

call_peaks() {
  local bam="$1"
  activate "$env_chipseq"
  macs2 callpeak -t "$bam" -g "$genome" --nomodel \
    --shift "$macs_shift" --extsize "$macs_extsize" \
    -n "$(basename "$bam" "_RD.nk.bam")" -B -f AUTO -q "$peak_q" \
    --outdir "$peak_dir" --SPMR
}

mkdir -p "$outdir" "$bw_dir" "$peak_dir"

echo "--- Mapping + dedup ---"
while read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  arr=($line)
  [ ${#arr[@]} -lt 3 ] && { echo "Error: need <sample> <fq1> <fq2>: $line"; exit 1; }
  sample="${arr[0]}"; fq1="${arr[1]}"; fq2="${arr[2]}"
  [ -f "$fq1" ] || { echo "Error: not found: $fq1"; exit 1; }
  [ -f "$fq2" ] || { echo "Error: not found: $fq2"; exit 1; }
  echo "Processing: $sample"
  map_track_a "$sample" "$fq1" "$fq2"
  map_track_b "$sample" "$fq1" "$fq2"
  to_bigwig "$sample"
done < "$sample_list"

echo "--- Peak calling ---"
shopt -s nullglob
bams=("${outdir}"/*_RD.nk.bam)
shopt -u nullglob
for bam in "${bams[@]}"; do
  call_peaks "$bam"
done

echo "Done"
