#!/usr/bin/env bash
# Time, peak memory and storage of FastQC-style statistics on simulated read
# pairs: seqqc over a morloc stream file and over FASTQ text, and
# FastQC 0.12.1 (the version nf-core/rnaseq 3.27.0 pins) over the FASTQ.
#
#   GENOME_SIZE   genome size in bases (default 100000000)
#   COVERAGES     read coverages (default "3 9")
#
# Columns: tool, input, pairs, seconds, peak_rss_mib, peak_disk_mib,
# written_mib.
set -euo pipefail
cd "$(dirname "$0")"
L=../..
P=$L/bioprojects
W=work
RES=results.tsv
HIST=history.tsv
MEASURE="python3 $P/tools/measure.py"
export MEASURE_DIRS="$PWD/$W/out"
SIZE=${GENOME_SIZE:-100000000}
CORES=$(getconf _NPROCESSORS_ONLN)
for x in "$P/synth/synth" "$P/trimmer/trimmer"; do
  [ -x "$x" ] || { echo "bench: $x missing" >&2; exit 1; }
done
[ -x ./bioseqqc-bench ] || morloc make -o bioseqqc-bench bench.loc
command -v fastqc >/dev/null || { echo "bench: fastqc is not on the PATH (bench/package.yaml)" >&2; exit 1; }
mkdir -p "$W/out"
$P/synth/synth -f packet -z 3 -o "$W/g.mlc" genome -n 200 -s "$SIZE" --seed 3
row() { local IFS=$'\t'; echo "$*" >> "$RES"; }
printf 'tool\tinput\tpairs\tseconds\tpeak_rss_mib\tpeak_disk_mib\twritten_mib\n' > "$RES"
for c in ${COVERAGES:-3 9}; do
  $P/trimmer/trimmer simulate "$W/g.mlc" -c "$c" --r1="$W/r1.fq" --r2="$W/r2.fq" --no-stdout
  n=$(( $(wc -l < "$W/r1.fq") / 4 ))
  ./bioseqqc-bench ingestFastqPairs --compression-level 3 "$W/r1.fq" "$W/r2.fq" "$W/r.mlc"
  row bioseqqc morloc-zstd3 "$n" "$($MEASURE ./bioseqqc-bench qcMlc "$W/r.mlc")"
  row bioseqqc fastq "$n" "$($MEASURE ./bioseqqc-bench qcFastq "$W/r1.fq" "$W/r2.fq")"
  rm -rf "$W/out"/*
  row fastqc-0.12.1 fastq "$n" "$($MEASURE fastqc -t "$CORES" -o "$W/out" "$W/r1.fq" "$W/r2.fq")"
  rm -rf "$W/out"/* "$W/r.mlc" "$W"/*.fq
done
commit=$(git -C /work/morloc-workspace/compiler/morloc rev-parse --short HEAD 2>/dev/null || echo NA)
[ -f "$HIST" ] || printf 'date\tcompiler\tlabel\t%s\n' "$(head -n 1 "$RES")" > "$HIST"
tail -n +2 "$RES" | while IFS= read -r line; do
  printf '%s\t%s\t%s\t%s\n' "$(date +%F)" "$commit" "${BENCH_LABEL:-run}" "$line" >> "$HIST"
done
rm -rf "$W"
cat "$RES"
