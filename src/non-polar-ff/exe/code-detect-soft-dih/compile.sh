#!/bin/bash

icc -O3 -o gen_soft_list gen_soft_list.cpp ff.cpp ../../../../opt/nlopt-2.4.2/.libs/libnlopt.a 

cp gen_soft_list ../../../../bin/
chmod g+rx ../../../../bin/gen_soft_list
