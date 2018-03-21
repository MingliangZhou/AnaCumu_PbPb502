#!/bin/sh

for ((I1=0;I1<2;I1+=1))
do
root -b -l <<EOF
.L ../CombSys.cxx+
CombSys($I1)
.q
EOF
done

rm -rf ../*.d
rm -rf ../*.so
