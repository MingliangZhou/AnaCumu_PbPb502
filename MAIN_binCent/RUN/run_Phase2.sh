#!/bin/sh

for ((I1=0;I1<1;I1+=1))
do
for ((I2=0;I2<=1;I2+=1))
do

root -b -l <<EOF
gSystem->Load(Cumu_cxx.so)
gSystem->Load(Phase2_cxx.so)
Phase2($I1,$I2)
.q
EOF

done
done
