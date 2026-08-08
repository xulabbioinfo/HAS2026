#!/bin/bash
## CUT&Tag paired-end pipeline: mapping (genome + spike-in) -> fragments -> dedup
##   -> peak calling -> spike-in scaled bigwig
## Usage: $0 [sample.list]
##   sample.list: <sample> <fastq1> <fastq2> per line

set -euo pipefail

## INPUTS ======================================================================
sample_list="${1:-./sample.list}"

genome_index=""           # bowtie2 genome index prefix (full path)
spikein_index=""          # bowtie2 spike-in index prefix (full path)

threads_bowtie2=25
threads_view=20
threads_sort=10
threads_bamcoverage=10

min_insert=10
max_insert=700
min_mapq=2
fragment_max=1000
bin_len=500
spikein_scale=1000000

genome="hs"
peak_q=0.01

outdir="mapping"
peak_dir="peakCalling"
## =============================================================================

bw_dir="${outdir}/bw"

[ -f "$sample_list" ] || { echo "Error: sample list not found: $sample_list"; exit 1; }
[ -n "$genome_index" ] || { echo "Error: genome_index is empty"; exit 1; }
[ -n "$spikein_index" ] || { echo "Error: spikein_index is empty"; exit 1; }
ls "${genome_index}".*.bt2 >/dev/null 2>&1 || { echo "Error: genome index not found: $genome_index"; exit 1; }
ls "${spikein_index}".*.bt2 >/dev/null 2>&1 || { echo "Error: spike-in index not found: $spikein_index"; exit 1; }

for cmd in bowtie2 samtools bedtools picard macs2 bamCoverage bc; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Error: $cmd not found in PATH"; exit 1; }
done

align() { # $1: index, $2: sample, $3: fq1, $4: fq2, $5: sam, $6: log
  bowtie2 --end-to-end --very-sensitive --no-unal --no-mixed --no-discordant \
    --phred33 -I "$min_insert" -X "$max_insert" -p "$threads_bowtie2" \
    -x "$1" -1 "$3" -2 "$4" -S "$outdir/$5" &> "$outdir/$6"
}

process_sample() { # $1: sample, $2: fq1, $3: fq2
  local sample="$1" fq1="$2" fq2="$3"

  align "$genome_index"  "$sample" "$fq1" "$fq2" "${sample}_bowtie2.sam"        "${sample}_bowtie2.txt"
  align "$spikein_index" "$sample" "$fq1" "$fq2" "${sample}_spikeIn_bowtie2.sam" "${sample}_spikeIn_bowtie2.txt"

  local seq_depth
  seq_depth=$(samtools view -F 4 "$outdir/${sample}_spikeIn_bowtie2.sam" | wc -l)
  echo $((seq_depth / 2)) > "$outdir/${sample}_bowtie2_spikeIn.seqDepth"

  samtools view -F 4 "$outdir/${sample}_bowtie2.sam" \
    | awk -F'\t' '{d=$9; print (d<0) ? -d : d}' \
    | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' \
    > "$outdir/${sample}_fragmentLen.txt"
  samtools view -h -@ "$threads_view" -q "$min_mapq" "$outdir/${sample}_bowtie2.sam" \
    > "$outdir/${sample}_filterminQualityScore.sam"

  samtools view -h -bS -F 4 "$outdir/${sample}_filterminQualityScore.sam" \
    > "$outdir/${sample}_bowtie2.mapped.bam"
  bedtools bamtobed -i "$outdir/${sample}_bowtie2.mapped.bam" -bedpe \
    > "$outdir/${sample}_bowtie2.bed"
  awk -v m="$fragment_max" '$1==$4 && $6-$2 < m' "$outdir/${sample}_bowtie2.bed" \
    > "$outdir/${sample}_bowtie2.clean.bed"
  cut -f 1,2,6 "$outdir/${sample}_bowtie2.clean.bed" \
    | sort -k1,1 -k2,2n -k3,3n > "$outdir/${sample}_bowtie2.fragments.bed"
  awk -v w="$bin_len" '{print $1, int(($2 + $3)/(2*w))*w + w/2}' \
    "$outdir/${sample}_bowtie2.fragments.bed" \
    | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' \
    | sort -k1,1V -k2,2n > "$outdir/${sample}_bowtie2.fragmentsCount.bin${bin_len}.bed"

  samtools sort "$outdir/${sample}_bowtie2.mapped.bam" -@ "$threads_sort" \
    -o "$outdir/${sample}_bowtie2.mapped.sort.bam"
  picard MarkDuplicates --INPUT "$outdir/${sample}_bowtie2.mapped.sort.bam" \
    --OUTPUT "$outdir/${sample}_bowtie2.mapped.sort_rmdup.bam" \
    --METRICS_FILE "$outdir/${sample}_bowtie2.mapped.sort_rmdup.log" \
    --CREATE_INDEX true --REMOVE_DUPLICATES true

  rm -f "$outdir/${sample}_bowtie2.sam" "$outdir/${sample}_spikeIn_bowtie2.sam" \
    "$outdir/${sample}_filterminQualityScore.sam"
}

call_peaks() {
  local sample="$1"
  macs2 callpeak -t "$outdir/${sample}_bowtie2.mapped.sort_rmdup.bam" \
    -g "$genome" -f BAMPE -n "${sample}_peak_q${peak_q}" --outdir "$peak_dir" \
    -q "$peak_q" -B --SPMR --keep-dup 1 2> "$peak_dir/${sample}_macs2Peak_summary.txt"
}

make_bigwig() {
  local sample="$1" seq_depth scale_factor
  seq_depth=$(cat "$outdir/${sample}_bowtie2_spikeIn.seqDepth")
  scale_factor=$(echo "$spikein_scale / $seq_depth" | bc -l)
  echo "Scaling factor for $sample: $scale_factor"
  bamCoverage --bam "$outdir/${sample}_bowtie2.mapped.sort_rmdup.bam" \
    -o "$bw_dir/${sample}_bowtie2.mapped.sort_rmdup.bw" \
    --binSize 1 --normalizeUsing RPKM --numberOfProcessors "$threads_bamcoverage" \
    --scaleFactor "$scale_factor"
}

mkdir -p "$outdir" "$peak_dir" "$bw_dir"

samples=()

echo "--- Mapping + processing ---"
while read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  arr=($line)
  [ ${#arr[@]} -lt 3 ] && { echo "Error: need <sample> <fq1> <fq2>: $line"; exit 1; }
  sample="${arr[0]}"; fq1="${arr[1]}"; fq2="${arr[2]}"
  [ -f "$fq1" ] || { echo "Error: not found: $fq1"; exit 1; }
  [ -f "$fq2" ] || { echo "Error: not found: $fq2"; exit 1; }
  echo "Processing: $sample"
  process_sample "$sample" "$fq1" "$fq2"
  samples+=("$sample")
done < "$sample_list"

[ ${#samples[@]} -gt 0 ] || { echo "Error: no samples in $sample_list"; exit 1; }

echo "--- Peak calling ---"
for sample in "${samples[@]}"; do
  call_peaks "$sample"
done

echo "--- Bigwig ---"
for sample in "${samples[@]}"; do
  make_bigwig "$sample"
done

echo "Done"
