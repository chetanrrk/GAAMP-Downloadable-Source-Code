#!/bin/bash


g++ -O3 -o check_lj modify_lj.cpp 

cp check_lj ../../../../bin/
chmod g+rx ../../../../bin/check_lj

