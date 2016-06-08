#!/bin/sh

#$ -l virtual_free=4G,mem_free=4G

echo started... $RAKE_TARGET
rake $RAKE_TARGET
echo ...finished.


