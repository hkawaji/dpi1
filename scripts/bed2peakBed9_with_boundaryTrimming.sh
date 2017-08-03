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
noise_subtraction_ratio=0

while getopts b:c:r: opt
do
  case ${opt} in
  b) bw=${OPTARG};;
  c) color=${OPTARG};;
  r) noise_subtraction_ratio=${OPTARG};;
  *) usage;;
  esac
done


if [ "${bw}" = "" ]; then usage; fi
if [ "${color}" = "" ]; then usage; fi
if [ "${noise_subtraction_ratio}" = "" ]; then usage; fi

awk --assign bw=${bw} --assign color=${color} --assign noise_subtraction_ratio=${noise_subtraction_ratio} \
'BEGIN{OFS="\t"}{

  chrom = $1; start = $2; stop = $3; name = $4; score = $5; strand = $6;

  ###
  ### for representative position
  ###
  rep_start = 0; rep_stop = 0; rep_max = 0;
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout",
  chrom, start, stop, bw )
  while ((command | getline ) > 0)
  {
    #total = total + ($4 * ($3 - $2))
    if ($4 > rep_max){
      rep_start = $2
      rep_stop  = $3
      rep_max   = $4
    }
  }
  close(command)

  ###
  ### for boundary trimming by noise subtraction
  ###
  boundary_start = 0; boundary_stop = 0 ; 
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout | sort -k1,1 -k2,2n ",
  chrom, start, stop, bw )
  while ((command | getline ) > 0)
  {
    if ( ( $4 - (rep_max * noise_subtraction_ratio) ) > 0  ){
      boundary_start = $2
      break
    }
  }
  close(command)
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout | sort -k1,1 -k2,2nr ",
  chrom, start, stop, bw )
  while ((command | getline ) > 0)
  {
    if ( ( $4 - (rep_max * noise_subtraction_ratio) ) > 0  ){
      boundary_stop = $3
      break
    }
  } 
  close(command)

  ###   
  ### get total counts on the new boundary
  ###   
  total = 0
  command = sprintf("bigWigToBedGraph -chrom=%s -start=%s -end=%s %s /dev/stdout",
  chrom, boundary_start, boundary_stop, bw )
  while ((command | getline ) > 0)
  {
    total = total + ($4 * ($3 - $2))
  }
  close(command)

  # debug
  #if (start != boundary_start) {print "xxxxx"}
  #if (stop  != boundary_stop)  {print "yyyy"}
  #if (score != total)  {print "zzz"}
  ###   
  ### print
  ###   
  if (total > 0)
  {
    # added this line at 13th Jan, 2012
    name=chrom":"boundary_start".."boundary_stop","strand

    print chrom, boundary_start, boundary_stop, name, total, strand,
          rep_start, rep_stop, color
  }
}'

