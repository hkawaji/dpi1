DPI1
====

Decomposition-based peak identification (DPI) - it finds peaks across a large number of TSS (transcription starting site) profiles. This set of scripts is tailored to produce CAGE peaks in FANTOM5 phase 1 (Forrest et. al., 2014, and its related papers), consisting of ~1000 TSS profiles. Note that this set of scripts is a typical patchwork, which depends on many programs and is tested only a few environments. I believe the FANTOM5 peaks are reasonably comprehensive for general analyses in human and mouse, and I won’t recommend to run this set of scripts in every projects. However, still, it would make a lot of sense to run this set of scripts for special purposes or cases. 


## Requirements 

  - ruby (https://www.ruby-lang.org)
  - R (http://cran.r-project.org/)
  - fastICA package (http://cran.r-project.org/web/packages/fastICA/index.html)
  - command line tools to operate bigWig files in http://hgdownload.cse.ucsc.edu/admin/
  - bedTools (https://code.google.com/p/bedtools/)
  - Unix/Linux with Grid Engine (developed with UGE)

## Installation

    % git clone https://github.com/hkawaji/dpi1.git

## How to run in this script

Step1.  set parameters in Rakefile

At least, these parameters in the Rakefile ( ${dpi_root}/Rakefile ) have to be set properly for your environment:

    % cd ${dpi_root}
    % head -10 Rakefile
    ...
    genome = "~/BEDTools/genomes/human.hg19.genome"
    ctss_path = "../in/ctss/*.ctss.bed.gz"
    out_path = "../dpiout/"
    ...

Step2. run

The shell script wrapper below run 'rake' consecutively via Grid Engine, which was required to handle thousands of CAGE data on the human genome. You could run individual commands without Grid Engine for a smaller set of data.

    % cd ${dpi_root}
    ./00run.sh >& 00run.err


## Input

  - CAGE read counts per CTSS (CAGE tag starting site - 5'end of CAGE reads) in BED format, gzipped.


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

2014 RIKEN, Japan. 

## Reference
A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. (accepted)


