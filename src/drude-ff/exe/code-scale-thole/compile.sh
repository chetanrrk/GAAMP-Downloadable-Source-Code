#!/bin/bash


g++ -O3 -o scale_thole scale-thole.cpp

cp scale_thole ../../../../bin
chmod g+rx ../../../../bin/scale_thole

