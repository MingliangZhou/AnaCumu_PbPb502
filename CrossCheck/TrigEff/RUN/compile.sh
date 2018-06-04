#!/bin/sh

root -b -l <<EOF
.L ../Event.cxx+
.L ../MultiCorr.cxx+
.L ../Phase1.cxx+
.q
EOF

rm -rf ../*.d
mv ../*.so .
