#!/bin/bash

g++ -O2 -o gen_xpsf gen_xpsf.cpp 

cp gen_xpsf ../../../../bin
chmod g+rx ../../../../bin/gen_xpsf
