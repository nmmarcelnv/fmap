#!/bin/bash

set -x

cat /proc/cpuinfo >test20.sh.cpuinfo
export OMP_NUM_THREADS=16
ln -sf ../dat/lattice.txt .
head -20 ../dat/3600.ang.dat >ang.dat.20
time ../src/fmap crowd.vdw 1.vdw 334 0.5970749 0.05 1.08 1.08 298 ang.dat.20 2.0 0.2 -9 1>test20.sh.out 2>test20.sh.err
grep -v "^  0.000000  0.000000  0.000000" test20.sh.out >fmap.out
time ../src/fmapdd fmap.out 0.05 298 2.0 0.2  >fmapdd.out
