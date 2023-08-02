#!/bin/bash

# compile
cd ${SRC_DIR}/src
make

# install
mkdir -p $PREFIX/bin
cp ${SRC_DIR}/src/hmmcnc $PREFIX/bin
cp ${SRC_DIR}/src/samToBed $PREFIX/bin
cp ${SRC_DIR}/src/viterbi $PREFIX/bin

