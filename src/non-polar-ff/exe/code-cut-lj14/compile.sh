#!/bin/bash

g++ -O3 -o cut_LJ14 cut-LJ14.cpp 

cp cut_LJ14 ../../../../bin/
chmod g+rx ../../../../bin/cut_LJ14
