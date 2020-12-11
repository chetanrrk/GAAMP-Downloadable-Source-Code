#!/bin/bash

g++ -O2 -o prep_cgenff extract_param.cpp ff.cpp 

cp prep_cgenff ../../../../bin/
chmod g+rx ../../../../bin/prep_cgenff
