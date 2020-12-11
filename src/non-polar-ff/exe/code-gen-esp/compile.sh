#!/bin/bash

g++ -O3 -o gen-esp gen-esp.cpp 

cp gen-esp ../../../../bin
chmod g+rx ../../../../bin/gen-esp
