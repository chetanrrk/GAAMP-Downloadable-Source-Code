#!/bin/bash


g++ -O3 -o scale_alpha scale-alpha.cpp

cp scale_alpha ../../../../bin
chmod g+rx ../../../../bin/scale_alpha

