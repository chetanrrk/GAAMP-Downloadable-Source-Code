#!/bin/bash

mpicxx -O3 -o qm-1d-scan QM_1D_scan.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp qm-1d-scan ../../../../bin/
chmod g+rx ../../../../bin/qm-1d-scan
