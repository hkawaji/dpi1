#!/usr/bin/env bash

function usage()
{
  cat <<EOF
usage: cat input.bed | $0 -b bigWig -c color

EOF
  exit 1;
}

bw=
color=

while getopts b:c: opt
do
  case ${opt} in
  b) bw=${OPTARG};;
  c) color=${OPTARG};;
  *) usage;;
  esac
done

if [ "${bw}" = "" ]; then usage; fi
if [ "${color}" = "" ]; then usage; fi

awk --assign bw=${bw} --assign color=${color} 'BEGIN{OFS="\t"}{
  chrom = $1
  start = $2
  stop  = $3
  name  = $4
  score  = $5
  strand = $6

  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout",
  chrom, start, stop, bw )

  rep_start = 0
  rep_stop  = 0
  rep_max   = 0
  total = 0
  while ((command | getline ) > 0)
  {
    total = total + ($4 * ($3 - $2))
    if ($4 > rep_max){
      rep_start = $2
      rep_stop  = $3
      rep_max   = $4
    }
  }
  close(command)

  score = total
  if (score > 0)
  {
    print chrom, start, stop, name, score, strand,
          rep_start, rep_stop, color
  }
}'

