#!/bin/sh

root -b -l <<EOF
.L ../Tool.cxx+
.L ../Event.cxx+
.L ../MultiCorr.cxx+
.q
EOF
