DPI1
====

Decomposition-based peak identification (DPI), which find peaks across a large number of TSS (transcription starting site) profiles 


## How to run in this script

  % ./dpi_run.sh -r

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
