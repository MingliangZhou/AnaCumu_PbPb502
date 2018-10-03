#!/bin/sh

for ((I1=0;I1<5;I1+=1))
do
root -b -l <<EOF
.L ../Plot.cxx+
Plot($I1)
.q
EOF
done

rm -rf ../*.d
rm -rf ../*.so
