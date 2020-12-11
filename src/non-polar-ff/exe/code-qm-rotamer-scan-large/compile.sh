#!/bin/bash

mpicxx -O3 -openmp -o qm-rotamer-scan-large QM_rotamer_scan.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp qm-rotamer-scan-large ../../../../bin/
chmod g+rx ../../../../bin/qm-rotamer-scan-large
