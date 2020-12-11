#!/bin/bash

g++ -O3 -o drude-xpsf gen_xpsf.cpp 

cp drude-xpsf ../../../../bin
chmod g+rx ../../../../bin/drude-xpsf
