#!/bin/bash


mpicxx -O2 -o org-rotamer-E 1D-rotamer-fitting.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp org-rotamer-E ../../../../bin
chmod g+rx ../../../../bin/org-rotamer-E

