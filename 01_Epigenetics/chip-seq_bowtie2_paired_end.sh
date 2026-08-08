#!/bin/bash
## ChIP-seq paired-end pipeline: mapping -> dedup -> bamqc -> peak calling -> bigwig
## Usage: $0 <sample.list> <info.list>
##   sample.list: <sample> <fastq1> <fastq2>
##   info.list: <treatment> <input>

set -euo pipefail

if [ $# -ne 2 ]; then
  echo "Usage: $0 <sample.list> <info.list>"
  exit 1
fi

sample_list="$1"
info_list="$2"
index=~/database/human/hg38_yuan_chip/UCSC/index/bowtie2/hg38
chrom_size=~/database/human/hg38_yuan_chip/UCSC/hg38.chrom.sizes
mapping_dir="mapping"

[ -f "$sample_list" ] || { echo "Error: $sample_list not found"; exit 1; }
[ -f "$info_list" ] || { echo "Error: $info_list not found"; exit 1; }

mkdir -p "$mapping_dir"

## mapping
echo "--- Mapping ---"
while read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  arr=($line)
  [ ${#arr[@]} -lt 3 ] && { echo "Error: need <sample> <fq1> <fq2>: $line"; exit 1; }
  sample="${arr[0]}"; fq1="${arr[1]}"; fq2="${arr[2]}"
  [ -f "$fq1" ] || { echo "Error: $fq1 not found"; exit 1; }
  [ -f "$fq2" ] || { echo "Error: $fq2 not found"; exit 1; }
  echo "Mapping: $sample"
  bowtie2 -p 25 -x "$index" -1 "$fq1" -2 "$fq2" \
    | samtools sort -O bam -@ 10 -o "${mapping_dir}/${sample}_bowtie2.bam" -
done < "$sample_list"

cd "$mapping_dir"

## flagstat
for bam in *.bam; do
  samtools flagstat -@ 2 "$bam" > "$(basename "$bam" .bam).flagstat"
done

## unique reads + remove duplicates
for bam in *_bowtie2.bam; do
  prefix="$(basename "$bam" .bam)"
  samtools view -@ 2 -F 4 -q 1 -b "$bam" | samtools sort -O bam -@ 2 -o "${prefix}.unique.bam" -
  picard MarkDuplicates --INPUT "${prefix}.unique.bam" \
    --OUTPUT "${prefix}.unique.bam_rmdup.bam" \
    --METRICS_FILE "${prefix}.unique.bam_rmdup.log" \
    --CREATE_INDEX true --REMOVE_DUPLICATES true
done

## bamqc
mkdir -p bam_qc
for bam in *_bowtie2.bam.unique.bam_rmdup.bam; do
  sample="$(basename "$bam" _bowtie2.bam.unique.bam_rmdup.bam)"
  Rscript /home/admin/software/phantompeakqualtools-master/run_spp_nodups.R \
    -c="$bam" -savp="./bam_qc/${sample}.plot.pdf" -out="./bam_qc/${sample}.qual" \
    -p=10 1> "./bam_qc/${sample}.Rout"
done

cd bam_qc
for f in *.qual; do head -n 3 "$f" >> result.txt; done

awk -F'\t' '
  BEGIN{OFS="\t"}
  {
    gsub(/^[ \t]+|[ \t]+$/, "", $1)
    split($1, arr, " ")
    sample = arr[1]
    split($3, frag, ",")
    print sample, frag[1]
  }' result.txt > sample.fragment

awk '
  NR==FNR { frag[$1]=$2; next }
  {
    if ($1 in frag && $2 in frag) {
      max = (frag[$1] > frag[$2]) ? frag[$1] : frag[$2]
      print $0, max
    }
  }' sample.fragment "../../$info_list" > peak_calling.list

cd ..

## peak calling
echo "--- Peak calling ---"
mkdir -p ../peak_calling
cp bam_qc/peak_calling.list ../peak_calling/
cd ../peak_calling

while read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  arr=($line)
  name="${arr[0]}"; input="${arr[1]}"
  treatment_bam="../mapping/${name}_bowtie2.bam.unique.bam_rmdup.bam"
  control_bam="../mapping/${input}_bowtie2.bam.unique.bam_rmdup.bam"
  [ -f "$treatment_bam" ] || { echo "Error: $treatment_bam not found"; exit 1; }
  [ -f "$control_bam" ] || { echo "Error: $control_bam not found"; exit 1; }
  echo "Calling peaks: $name"
  macs2 callpeak -t "$treatment_bam" -c "$control_bam" \
    -f BAMPE -g hs -n "$name" -B --SPMR -q 0.01 \
    --outdir . 2> "${name}_macs2Peak_summary.txt"
done < peak_calling.list

## bigwig
echo "--- Bigwig ---"
mkdir -p bw && cd bw
while read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  arr=($line)
  name="${arr[0]}"
  macs2 bdgcmp -t "../${name}_treat_pileup.bdg" -c "../${name}_control_lambda.bdg" \
    -o "./${name}_subtract.bdg" -m subtract
  awk '{OFS="\t"; if ($4 < 0) {print $1, $2, $3, 0} else {print $0}}' "./${name}_subtract.bdg" \
    > "./${name}_subtract_nonegative.bdg"
  /home/admin/software/bdgcmp2bw.sh "${name}_subtract_nonegative.bdg" "$chrom_size"
done < ../peak_calling.list

echo "Done: $(date +"%Y-%m-%d %H:%M:%S")"
