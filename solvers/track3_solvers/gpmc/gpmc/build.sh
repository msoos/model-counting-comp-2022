#!/bin/sh

mkdir -p build
cd build

if [ $1 = "r" ]; then
  echo "Create Makefile for Release..."
  cmake -DCMAKE_BUILD_TYPE=Release .. \
  && make \
  && cp -p gpmc ../../bin/ 

elif [ $1 = "clean" ]; then
  make clean

fi

