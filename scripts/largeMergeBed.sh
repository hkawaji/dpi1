#!/bin/sh

distance=20
tmpdir=.

function usage()
{
  cat <<EOF
usage: cat input.bed | $0 -g genome [-d min_distance] [-t tmpdir]

Note: The input BED (-i) file must be grouped by chromosome.
      A simple "sort -k 1,1 <BED> > <BED>.sorted" will suffice.

EOF
  exit 1;
}

while getopts t:d:g: opt
do
  case ${opt} in
  t) tmpdir=${OPTARG};;
  d) distance=${OPTARG};;
  g) genome=${OPTARG};;
  *) usage;;
  esac
done

if [ "${genome}" = "" ]; then usage; fi


extend_bp=$(expr $distance / 2)

#sort -k1,1 --temporary-directory=${tmpdir} \

slopBed -b ${extend_bp} -i stdin -g ${genome} \
| genomeCoverageBed -i stdin -g ${genome} -bga \
| grep --perl-regexp "\t0$" \
| complementBed -i stdin -g ${genome} \
| slopBed -b -${extend_bp} -i stdin -g ${genome} \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2".."$3}' \

