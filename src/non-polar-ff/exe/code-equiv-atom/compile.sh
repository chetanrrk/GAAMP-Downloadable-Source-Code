#!/bin/bash

g++ -O2 -o equiv_atom equiv_atom.cpp

cp equiv_atom ../../../../bin
chmod g+rx ../../../../bin/equiv_atom
