#!/bin/bash

g++ -O3 -o check-b0-theta0 check_b0_theta0.cpp ff.cpp 

cp check-b0-theta0 ../../../../bin/
chmod g+rx ../../../../bin/check-b0-theta0
