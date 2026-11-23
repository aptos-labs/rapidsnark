#!/bin/bash

set -e

ROOTPATH=$(pwd)
INSTALL_PREFIX="${1:-$ROOTPATH/package}"
if [[ $# -gt 0 ]]; then
  shift
fi

cd depends/tbb/oneTBB
# Create binary directory for out-of-source build
rm -rf build && mkdir -p build && cd build
# Configure: customize CMAKE_INSTALL_PREFIX and disable TBB_TEST to avoid tests build
cmake -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" -DCMAKE_INSTALL_LIBDIR=lib -DTBB_TEST=OFF -DTBB_STRICT=OFF "$@" ..
# Build
cmake --build . --parallel
# Install
cmake --install .

cd ../../..
