#!/bin/sh

#$ -l virtual_free=4G,mem_free=4G

if [ "${PATH_OVERRIDE}" != "" ]; then
  PATH=${PATH_OVERRIDE}
fi

echo started... $RAKE_TARGET
rake $RAKE_TARGET
echo ...finished.


