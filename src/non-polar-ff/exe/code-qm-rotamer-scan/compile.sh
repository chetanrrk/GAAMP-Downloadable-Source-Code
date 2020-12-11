#!/bin/bash


mpicxx -O3 -o qm-rotamer-scan QM_rotamer_scan.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp qm-rotamer-scan ../../../../bin/
chmod g+rx ../../../../bin/qm-rotamer-scan
