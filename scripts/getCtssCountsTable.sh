#!/usr/bin/env bash

# setup
set -e
tmpdir=$(mktemp -d)
trap "[[ $tmpdir ]] && rm -rf $tmpdir" 0 1 2 3 15

function usage()
{
  cat <<EOF

$0 obtain CTSS values as a matrix, in which
rows and columns indicate the input file names

  usage: $0 -c chrom -s start -e end -i input_file_list

columns orderd from start to the end per each base pair.

EOF
  exit 1;
}

#chrom=chr9
#start=2015336
#stop=2015370
#input_file_list=list.txt

while getopts c:s:e:i:n:p: opt
do
  case ${opt} in
  c) chrom=${OPTARG};;
  s) start=${OPTARG};;
  e) stop=${OPTARG};;
  i) input_file_list=${OPTARG};;
  p) parallel_jobs=${OPTARG};;
  *) usage;;
  esac
done

if [ "${chrom}" = "" ]; then usage; fi
if [ "${start}" = "" ]; then usage; fi
if [ "${stop}" = "" ]; then usage; fi
if [ "${input_file_list}" = "" ]; then usage; fi

for infile in $( cat ${input_file_list} )
do
  fname=$(basename $infile)
  bigWigToBedGraph -chrom=${chrom} -start=${start} -end=${stop} ${infile} /dev/stdout \
  | awk --assign fname=$fname 'BEGIN{OFS="\t"}{for (i=$2;i<$3;i++){print fname, $1":"i".."i+1,$4}}'
done | R --slave -e  'library(tidyr);read.table(file("stdin"),as.is=T) %>% tidyr::spread("V1", "V3", fill=0, convert=FALSE) %>% write.table(quote=F,sep="\t",row.names=F)'

