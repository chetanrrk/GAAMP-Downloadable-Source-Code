#!/bin/bash

g++ -O2 -o exclude_H2 exclude-H2.cpp 

cp exclude_H2 ../../../../bin
chmod g+rx ../../../../bin/exclude_H2
