#!/bin/bash

mpicxx -O3 -o mm_pes_large mm_pes.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp mm_pes_large ../../../../bin/
chmod g+rx ../../../../bin/mm_pes_large

