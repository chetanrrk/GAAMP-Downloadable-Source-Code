#!/bin/bash

mpicxx -O2 -o qm-1d-scan-para QM_1D_scan.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp qm-1d-scan-para ../../../../bin/
chmod g+rx ../../../../bin/qm-1d-scan-para

