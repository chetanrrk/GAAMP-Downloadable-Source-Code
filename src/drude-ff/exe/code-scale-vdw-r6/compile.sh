#!/bin/bash


g++ -O3 -o scale_vdw_r6 scale-vdw-r6.cpp

cp scale_vdw_r6 ../../../../bin
chmod g+rx ../../../../bin/scale_vdw_r6

