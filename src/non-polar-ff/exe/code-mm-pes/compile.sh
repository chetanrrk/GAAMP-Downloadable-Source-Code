#!/bin/bash

mpicxx -O3 -o mm_pes mm_pes.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp mm_pes ../../../../bin
chmod g+rx ../../../../bin/mm_pes
