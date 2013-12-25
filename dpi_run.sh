#!/bin/sh

### param

target_tag_cluster_in_bed='example_tag_cluster.bed'
path_to_bigWig='in/'
pattern_of_bigWig=".bw$"
exclude_prefix='NA'
gaussian_window_size_half=5
n_comp_upper_bound=5
length_to_decompose=50
noise_subtraction_ratio=0.1

### R code

code=$(cat <<EOF
source("dpi_core.R");
  main(
    infile_base="${target_tag_cluster_in_bed}",
    path="${path_to_bigWig}" ,
    pattern='${pattern_of_bigWig}' ,
    exclude_prefix=${exclde_prefix} ,
    gaussian_window_size_half= ${gaussian_window_size_half},
    n.comp.upper_bound = ${n_comp_upper_bound},
    length_to_decompose = ${length_to_decompose} ,
    noise_subtraction_ratio = ${noise_subtraction_ratio}
  )
EOF
)

function usage()
{
  cat <<EOF

DPI1
====

Decomposition-based peak identification (DPI), which find peaks across a large number of TSS (transcription starting site) profiles 


## How to run in this script

    ./dpi_run.sh -r

'dpi_core.R' is the core code, and this shell script is a wrapper.

## Requirements 

  - R
  - fastICA package (http://cran.r-project.org/web/packages/fastICA/index.html)
  - bigWigToBedGraph in jksrc.zip (http://hgdownload.cse.ucsc.edu/admin/)
  - bedTools (https://code.google.com/p/bedtools/)

## Input

  - tag cluster to be segregated (infile_base)
  - CAGE read counts per CTSS (5'end of CAGE reads), formatted in bigWig 
    (*.fwd.bw and *.rev.bw for forward and reverse strands)

## Output

  - bed file


## Author

Hideya Kawaji


## Reference
A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. (submitted)


EOF
  exit 1;
}


while getopts hr opt
do
  case ${opt} in
  r) run=yes;;
  *) usage;;
  esac
done

if [ "${run}" != "yes" ]; then usage; fi


# run the code
echo $code \
| R --slave \
| grep -v '^#' \
| mergeBed -s -i stdin \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2".."$3","$4,0,$4}'


