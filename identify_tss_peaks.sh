#!/bin/bash

function usage()
{
  cat <<EOF

Decomposition-based peak identification (DPI)
=============================================

This tool identifies a set of reference peaks in in genome-wide TSS (transcription
starting site) profiles obtained from diverse range of biological states. In
particular, this is developed for the FANTOM5 project (Nature 507, 462-467, 2014
and its related papers) that produced CAGE profiles in more than one thousand of
biological states in human and mouse, and the analysis results have been
available at http://fantom.gsc.riken.jp/5/ as well as the article above.

In general, I recommend to use the FANTOM5 peak set for CAGE analysis, rather
than runnning this tool for your own set, so that you can compare your results
with others. I would recommend to run this to identify a set of TSS peaks in
other organisms than human and mouse, or very unique biological condition not
included in the FANTOM5 sample collecttion.


SYNOPSIS
--------
usage: $0 -g GENOME -i CTSS_FILES -o OUT_PATH [ -d DECOMPOSITION(y/N) ]


OPTIONS
-------
* -g GENOME: chrom_size file in BEDTools https://github.com/arq5x/bedtools

* -i CTSS_FILES: CAGE profiles as counts per CTSS (CAGE tag starging sites) in BED
  format. Examples can be found in the FANTOM5 web resource, such as
  http://fantom.gsc.riken.jp/5/datafiles/latest/basic/human.cell_line.LQhCAGE/.
  Specify them by "-i './human.cell_line.hCAGE/*bed.gz'"

* -d DECOMPOSITION(y/N): Specify if decomposition (based on ICA) has to be applied
  or not. The 'decomposition' step is incorporated to cope with a large number of (~1,000)
  diversed biological states such as FANTOM5 phase1 (Nature 507, 462-467, 2014). The
  step requires specific configuration of a computer cluster (with grid engine) to
  perform the calculation efficiently, however, it would not be required for usual
  cases using smaller or homogenous data set. As default, this script run without the
  decomposition step, which requires ruby, R, bigWig tools, and bedtools but not
  GridEngine nor R/fastICA. You can tun it on by setting "-d Y" to include the
  decomposition step after setting up your computer properly.

* -o OUT_PATH: all of the resulting files, including intermediate files. See 
  Nature 507, 462-467, 2014 for the definition of thresholding ("robust" and
  "permissive"). 'decompose_smoothing' is the identified peaks after decomposition,
  and 'spi_merged' is the one without decomposition. Following files would be the
  final results after thresholding:
	- outPooled/tc.[decompose_smoothing|spi_merged].bed.gz - all peaks defined.
	- outPooled/tc.[decompose_smoothing|spi_merged].ctssMaxCounts11_ctssMaxTpm1.bed.gz - 'robust' threshold
	- outPooled/tc.[decompose_smoothing|spi_merged].ctssMaxCounts3.bed.gz - 'permissive' threshold


AUTHOR
------
KAWAJI, Hideya


COPYRIGHT
---------
2014 RIKEN, Japan.


REFERENCE
---------
* A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. Nature 507, 462-467, 2014
* FANTOM5 web resource ( http://fantom.gsc.riken.jp/5/ )


EOF
  exit 1;
}


RAKE_OPT=""
decomposition=N
while getopts d:g:i:o:b:m:n:l:r:u:y opt
do
  case ${opt} in
  d) export decomposition=${OPTARG};;
  g) export dpi_genome=${OPTARG};;
  i) export dpi_ctss_path="${OPTARG}";;
  o) export dpi_out_path=${OPTARG};;
  b) export dpi_cluster_dist=${OPTARG};;
  m) export dpi_ctss_max_counts_threshold=${OPTARG};;
  n) export dpi_noise_subtraction_ratio=${OPTARG};;
  l) export dpi_length_to_decompose=${OPTARG};;
  r) export dpi_count_to_decompose=${OPTARG};;
  u) export dpi_n_comp_upper_bound=${OPTARG};;
  y) RAKE_OPT="${RAKE_OPT} --dry-run " ;;
  *) usage;;
  esac
done

if [ "${dpi_genome}" = "" ]; then usage; fi
if [ "${dpi_ctss_path}" = "" ]; then usage; fi
if [ "${dpi_out_path}" = "" ]; then usage; fi

export curr_dir=$(pwd)
installed_dir=$(dirname $0)

### prep for qsub
qsub_v="dummy=dummy"
for k in curr_dir installed_dir decomposition \
         dpi_genome dpi_ctss_path dpi_out_path \
         dpi_cluster_dist dpi_ctss_max_counts_threshold  \
         dpi_noise_subtraction_ratio dpi_length_to_decompose \
         dpi_count_to_decompose dpi_n_comp_upper_bound
do
  v='$'$k
  v=$(eval echo \"$v\")
  if [ "$v" != "" ]; then
    qsub_v="${k}=${v},${qsub_v}"
  fi
done
echo "qsub_v:$qsub_v" 1>&2

case ${decomposition} in
  "n"|"N")
    (
    cd $installed_dir
    rake ${RAKE_OPT} clean_all
    rake ${RAKE_OPT} prep
    rake ${RAKE_OPT} bw
    rake ${RAKE_OPT} tc
    rake ${RAKE_OPT} tc_short
    rake ${RAKE_OPT} tc_long
    rake ${RAKE_OPT} spi
    )
  ;;
  "y"|"Y")
    (
    cd $installed_dir
    ### prep ###
    rake ${RAKE_OPT} clean_all
    rake ${RAKE_OPT} prep

    ### 01 ###
    for X in $(rake --dry-run bw 2>&1 | grep Execute | cut -f 5 -d ' ' | grep -v '^bw$')
    do 
      qsub -o $curr_dir -e $curr_dir -v "${qsub_v}" -N xqsub$$_01_bw ./scripts/qsub.sh $X
    done

    ### 02 - 05 ###
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_02_pool              -hold_jid xqsub$$_01_bw   ./scripts/qsub.sh  pool_ctss
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_03_tag_clustering    -hold_jid xqsub$$_02_pool ./scripts/qsub.sh  tc
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_04_tc_short          -hold_jid xqsub$$_03_tag_clustering ./scripts/qsub.sh tc_short
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_05_tc_long -sync y   -hold_jid xqsub$$_03_tag_clustering ./scripts/qsub.sh tc_long

    ### 06 ###
    for X in $(rake --dry-run decompose_smoothing 2>&1 | grep Execute | cut -f 5 -d ' ' | grep  '.bed$' )
    do
      qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_06_dpi_each -hold_jid xqsub$$_05_tc_long ./scripts/qsub.sh $X
    done

    ### 07 - 09 ###
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_07_comp_ctss         -hold_jid xqsub$$_06_dpi_each ./scripts/qsub.sh comp_ctss
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_08_final             -hold_jid xqsub$$_06_dpi_each,xqsub$$_04_tc_short ./scripts/qsub.sh final
    qsub -o $curr_dir -e $curr_dir -v ${qsub_v} -N xqsub$$_09_threshold -sync y -hold_jid xqsub$$_08_final ./scripts/qsub.sh threshold
    )
  ;;
  *)
    usage;;
esac


