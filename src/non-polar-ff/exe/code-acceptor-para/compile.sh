#!/bin/bash


icc -O3 -o acceptor_para mol_H_acceptor.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp acceptor_para ../../../../bin
chmod g+rx ../../../../bin/acceptor_para
