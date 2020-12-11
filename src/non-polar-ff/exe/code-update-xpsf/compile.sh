#!/bin/bash

g++ -O3 -o update-xpsf update-xpsf.cpp 

cp update-xpsf ../../../../bin/
chmod g+rx ../../../../bin/update-xpsf
