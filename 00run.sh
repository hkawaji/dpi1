#!/bin/sh

### prep ###
rake prep

### 01 ###
for X in $(rake --dry-run bw 2>&1 | grep Execute | cut -f 5 -d ' ' | grep -v '^bw$')
do
  qsub -N xqsub$$_01_bw ./scripts/qsub.sh $X
done

### 02 - 05 ###
qsub -N xqsub$$_02_pool              -hold_jid xqsub$$_01_bw   ./scripts/qsub.sh  pool_ctss
qsub -N xqsub$$_03_tag_clustering    -hold_jid xqsub$$_02_pool ./scripts/qsub.sh  tc
qsub -N xqsub$$_04_tc_short          -hold_jid xqsub$$_03_tag_clustering ./scripts/qsub.sh tc_short
qsub -N xqsub$$_05_tc_long -sync y   -hold_jid xqsub$$_03_tag_clustering ./scripts/qsub.sh tc_long

### 06 ###
#for X in $( ls ../out_pooled/tc.long/* | grep -v .bed$ )
for X in $(rake --dry-run decompose_smoothing 2>&1 | grep Execute | cut -f 5 -d ' ' | grep  '.bed$' )
do
  qsub -N xqsub$$_06_dpi_each -hold_jid xqsub$$_05_tc_long ./scripts/qsub.sh $X
done

### 07 - 09 ###
qsub -N xqsub$$_07_comp_ctss         -hold_jid xqsub$$_06_dpi_each ./scripts/qsub.sh comp_ctss
qsub -N xqsub$$_08_final             -hold_jid xqsub$$_06_dpi_each,xqsub$$_04_tc_short ./scripts/qsub.sh final
qsub -N xqsub$$_09_threshold -sync y -hold_jid xqsub$$_08_final ./scripts/qsub.sh threshold
