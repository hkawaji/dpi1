Decomposition-based peak identification (DPI)
=============================================

This tool identifies a set of reference peaks in in genome-wide TSS (transcription
starting site) profiles obtained from diverse range of biological states. In
particular, this is developed for the FANTOM5 project (Nature 507, 462-467, 2014
and its related papers) that produced CAGE profiles in more than one thousand of
biological states in human and mouse, and the analysis results have been
available at http://fantom.gsc.riken.jp/5/ as well as the article above.

In general, I recommend to use the FANTOM5 peak set for CAGE analysis, rather
than runnning this tool for your own data, so that you can compare your results
with others. I would recommend to run this tool to identify a set of TSS peaks in
other organisms than human and mouse, or very unique biological conditions not
included in the FANTOM5 sample collecttion.


Requirements 
------------

  - bash (https://www.gnu.org/software/bash/)
  - ruby (https://www.ruby-lang.org)
  - R (http://cran.r-project.org/)
  - R/fastICA package (http://cran.r-project.org/web/packages/fastICA/)
  - R/tidyr package (https://cran.r-project.org/web/packages/tidyr/)
  - command line tools to operate bigWig files in http://hgdownload.cse.ucsc.edu/admin/
  - bedTools (https://code.google.com/p/bedtools/)
  - Unix/Linux with Grid Engine (developed with UGE)

Installation
------------

    % git clone https://github.com/hkawaji/dpi1.git

How to run
-----------

    % cd ${INSTALLED_DIR}
    % ./identify_tss_peaks.sh [options]

Follow the usage described in the message.


Update (branches)
-----------------
* beta3 (current master) (Aug, 2018)
  - updated for platform compatibility
* beta2 (June, 2016)
  - updated for compatibility to the recent version of bedtools (v2.25.0) and jksrc (v333).
  - updated for portability of user environment
* beta (Sep, 2014)
  - a wrapper script is updated as 'identify_tss_peaks.sh',
    enabling to run the steps in any directory.
  - set 'non-decomposition' mode as default
* alpha (April, 2014)
  - a wrapper script '00run.sh' perform all the requried steps.
* core (Dec, 2013):
  - only the core script of R


Author
------
Hideya Kawaji


Copyright
---------
2014-2017 RIKEN, Japan. 


Reference
---------
* A promoter level mammalian expression atlas, Forrest A, Kawaji H, Rehli M, et al. Nature 507, 462-467, 2014


