#!/bin/sh

for ((I1=0;I1<2;I1+=1))
do
for ((I2=0;I2<=5;I2+=1))
do
root -b -l <<EOF
.L ../Cumu.cxx+
.L ../Phase2.cxx+
Phase2($I1,$I2)
.q
EOF
done
done

rm -rf ../*.d
rm -rf ../*.so
