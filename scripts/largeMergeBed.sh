#!/usr/bin/env bash

distance=20
tmpdir=.

function usage()
{
  cat <<EOF
usage: cat input.bed | $0 [-d min_distance] 

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

###
### The steps has been largely simplified, owing the speed-up
### of mergeBed in bedtools.
###
### several parameters remains for backward compatibility,
### even if they are not required now.
###

sort -k1,1 -k2,2n \
| mergeBed -d $distance \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2".."$3}' \


