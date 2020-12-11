#!/bin/bash


g++ -O3 -o drude-ff modify_rtf.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp drude-ff ../../../../bin
chmod g+rx ../../../../bin/drude-ff
