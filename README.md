dpi1
====

decomposition-based peak identification, which find peaks across a large number of TSS (transcription starting site) profiles 

This is a sample script to perform decomposition-based peak identification (DPI).
The core code is dpi_core.R.


## how to run in this script

  % ./dpi_run.sh -r

## requirements 

  - R
  - fastICA package (http://cran.r-project.org/web/packages/fastICA/index.html)
  - bigWigToBedGraph in jksrc.zip (http://hgdownload.cse.ucsc.edu/admin/)
  - bedTools (https://code.google.com/p/bedtools/)

## input

  - tag cluster to be segregated (infile_base)
  - CAGE read counts per CTSS (5'end of CAGE reads), formatted in bigWig 
    (*.fwd.bw and *.rev.bw for forward and reverse strands)

## output

  - bed file


## Author

Hideya Kawaji


## reference
A promoter level mammalian expression atlas
Forrest A, Kawaji H, Rehli M, et al. (submitted)
