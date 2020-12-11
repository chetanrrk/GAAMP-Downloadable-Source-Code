#!/bin/bash


icc -O3 -o acceptor mol_H_acceptor.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp acceptor ../../../../bin/
chmod g+rx ../../../../bin/acceptor
