#!/bin/sh

#$ -cwd
#$ -l virtual_free=4G,mem_free=4G

source ~/.bashrc
echo $1 started...
rake $1
echo ...finished.


