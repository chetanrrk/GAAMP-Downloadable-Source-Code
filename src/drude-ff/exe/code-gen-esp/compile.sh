#!/bin/bash

g++ -O3 -o gen-esp-drude gen-esp.cpp 

cp gen-esp-drude ../../../../bin
chmod g+rx ../../../../bin/gen-esp-drude

