#!/bin/bash

ifort -O2 -o cgrid ConnollyGrid.f Surface.f -static

cp cgrid ../../../../bin/
chmod g+rx ../../../../bin/cgrid
