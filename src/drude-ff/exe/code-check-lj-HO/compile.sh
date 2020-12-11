#!/bin/bash


g++ -O3 -o check_lj-drude modify_lj.cpp 

cp check_lj-drude ../../../../bin
chmod g+rx ../../../../bin/check_lj-drude

