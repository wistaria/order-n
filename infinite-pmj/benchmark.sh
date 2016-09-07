#!/bin/sh

SIZES="4 8 16 32 64 128 256 512 1024 2048 4096"

SEED=12345
P_FERRO=0.5
T=1.0
MCS=1024
EPSILON=0.1

for N in ${SIZES}; do
  SPEED_NAIVE=$(./naive_mc ${SEED} ${N} ${P_FERRO} ${T} ${MCS} 2>&1 | grep Speed | awk '{print $3}')
  SPEED_ON=$(./order_n ${SEED} ${N} ${P_FERRO} ${T} ${MCS} ${EPSILON} 2>&1 | grep Speed | awk '{print $3}')
  echo ${N} ${SPEED_NAIVE} ${SPEED_ON}
done
