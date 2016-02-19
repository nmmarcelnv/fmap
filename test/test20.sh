#!/bin/bash

set -x

export MV2_USE_AFFINITY=0
export MV2_ENABLE_AFFINITY=0
export VIADEV_USE_AFFINITY=0
export VIADEV_ENABLE_AFFINITY=0

tail -25 /proc/cpuinfo >test20.sh.cpuinfo
free >>test20.sh.cpuinfo

export OMP_NUM_THREADS=16
export MIC_OMP_NUM_THREADS=240
ln -sf ../dat/lattice.txt .
head -20 ../dat/3600.ang.dat >ang.dat.20
time ../src/fmap crowd.vdw 1.vdw 334 0.5970749 0.05 1.08 1.08 298 ang.dat.20 2.0 0.2 -9 1>test20.sh.out 2>test20.sh.err
grep -e Clck -e Time test20.sh.err
grep -v "^  0.000000  0.000000  0.000000" test20.sh.out >fmap.out
time ../src/fmapdd fmap.out 0.05 298 2.0 0.2  >fmapdd.out
time ../src/fmapdd.mic fmap.out 0.05 298 2.0 0.2  >fmapdd.mic.out
grep v+e test20.sh.err >test20.sh.err.v+e
~/bin/colstat.sh test20.sh.err.v+e 3|awk '{printf("%16e\n",$1)}'
echo 1.0 1.0| ../src/matevscl mat.bin 4001 0.02 -40.0 298 
