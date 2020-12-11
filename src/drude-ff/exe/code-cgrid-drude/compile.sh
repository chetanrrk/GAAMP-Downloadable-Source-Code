#!/bin/bash

ifort -O2 -o drude-cgrid ConnollyGrid.f Surface.f -static

cp drude-cgrid ../../../../bin/
chmod g+rx ../../../../bin/drude-cgrid
