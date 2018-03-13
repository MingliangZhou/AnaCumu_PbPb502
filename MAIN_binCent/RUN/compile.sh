#!/bin/sh

root -b -l <<EOF
.L ../Tool.cxx+
.L ../Event.cxx+
.L ../MultiCorr.cxx+
.L ../Phase1.cxx+
.L ../Cumu.cxx+
.L ../Phase2.cxx+
.q
EOF

rm -rf ../*.d
mv ../*.so .
