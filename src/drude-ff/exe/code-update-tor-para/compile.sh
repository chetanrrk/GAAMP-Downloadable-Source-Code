#!/bin/bash


icc -O3 -o update-tor-para-drude update-torsion-para.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp update-tor-para-drude ../../../../bin/
chmod g+rx ../../../../bin/update-tor-para-drude
