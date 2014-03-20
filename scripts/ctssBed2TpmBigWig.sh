#!/bin/sh


### prep
function usage()
{
  cat <<EOF
usage: $0 -i INFILE.ctss.bed -o OUT_PREFIX -g GENOME
EOF
  exit 1;
}

while getopts i:o:g: opt
do
  case ${opt} in
  i) infile=${OPTARG};;
  o) outprefix=${OPTARG};;
  g) genome=${OPTARG};;
  *) usage;;
  esac
done

if [ "${infile}" = "" ]; then usage; fi
if [ "${outprefix}" = "" ]; then usage; fi
if [ "${genome}" = "" ]; then usage; fi




fwd=${outprefix}.fwd.bw
rev=${outprefix}.rev.bw
tmpfile=${outprefix}.tmp.bg

sum=$( gunzip -c ${infile} | awk 'BEGIN{sum=0}{sum=sum+$5}END{print sum}' )

### fwd
gunzip -c ${infile} \
| grep ^chr \
| grep +$ \
| awk --assign sum=${sum} 'BEGIN{OFS="\t";CONVFMT="%f"}{print $1,$2,$3, 1e6 * $5 / sum}' \
| sort -k1,1 -k2,2n \
> ${tmpfile} 

bedGraphToBigWig ${tmpfile} ${genome} ${fwd}

### rev
gunzip -c ${infile} \
| grep ^chr \
| grep -v +$ \
| awk --assign sum=${sum} 'BEGIN{OFS="\t";CONVFMT="%f"}{print $1,$2,$3, 1e6 * $5 / sum}' \
| sort -k1,1 -k2,2n \
> ${tmpfile} 
bedGraphToBigWig ${tmpfile} ${genome} ${rev}

rm -f ${tmpfile}

