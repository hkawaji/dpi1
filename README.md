DPI1
====

Decomposition-based peak identification (DPI), which find peaks across a large number of TSS (transcription starting site) profiles. This set of scripts is the one used in a paper (Forrest et. al., 2014) and its related papers.

## How to run in this script

Step1.  preparation - set parameters in Rakefile

At least, these parameters have to be set properly:

    genome = "~/BEDTools/genomes/human.hg19.genome"
    ctss_path = "../in/ctss/*.ctss.bed.gz"
    out_path = "../dpiout/"

Step2. run

The shell script wrapper below run 'rake' consecutively via Grid Engine, which was required to handle thousands of CAGE data on the human genome. You could run individual commands without Grid Engine for a smaller set of data.

    ./00run.sh >& 00run.err


## Requirements 

  - ruby (https://www.ruby-lang.org)
  - R (http://cran.r-project.org/)
  - fastICA package (http://cran.r-project.org/web/packages/fastICA/index.html)
  - command line tools to operate bigWig files in http://hgdownload.cse.ucsc.edu/admin/
  - bedTools (https://code.google.com/p/bedtools/)
  - Grid Engine (developed with UGE)

## Input

  - CAGE read counts per CTSS (5'end of CAGE reads) in BED format, gzipped.

## Output

Many files are generated under ${out_path}. Main results are:

  - ${out_path}/tc.decompose_smoothing_merged.ctssMaxCounts11_ctssMaxTpm1.bed.gz -  robust peaks
  - ${out_path}/tc.decompose_smoothing_merged.ctssMaxCounts3.bed.gz - permissive peaks

Additional useful files:

  - ${out_path}/tc.decompose_smoothing_merged.bed.gz - all peaks
  - ${out_path}/outCounts - bigWig files of CTSS counts
  - ${out_path}/outTpm - bigWig files of CTSS activity (normalized counts)


## Author

Hideya Kawaji


## Copyright

2014 RIKEN, Japan. ALL RIGHTS RESERVED. 


## Reference
A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. (accepted)


