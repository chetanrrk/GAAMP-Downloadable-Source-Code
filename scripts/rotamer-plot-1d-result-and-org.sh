#!/bin/bash

export GDFONTPATH=/home/huanglei/tools/gnuplot-install/fonts

echo "" > result-1D.html

NINC=1
NumTorsion=1
for (( ; ; ))
do
  if [ -e ./1d-qm-mm-$NumTorsion.dat ]; then
    title=`head -n $NumTorsion soft-dih-list.txt | tail -n 1 | awk '{print "Dihedral  " $1,"-",$2,"-",$3,"-",$4}'`
    echo "set terminal png medium font \"arial,26\" size 1200,800" > plot.sh
    echo "set title \"$title\"" >> plot.sh
    echo "set xlabel \"phi-$NumTorsion\" font \"arial,26\"" >> plot.sh
    echo "set ylabel \"E\" font \"arial,26\"" >> plot.sh
    echo "set output \"plot-1d-$NumTorsion.png\"" >> plot.sh
    echo "plot \"1d-qm-mm-$NumTorsion.dat\" using 1:2 t \"QM\" with linespoints pt 6 lw 2 ps 2 lc rgb \"black\", \"1d-qm-mm-$NumTorsion.dat\" using 1:3 t \"Fitted\" with linespoints pt 1 lw 2 ps 2 lc 1, \"org-1d-$NumTorsion.dat\" using 1:3 t \"org\" with linespoints pt 2 lw 2 ps 2 lc 2" >> plot.sh
  
    /home/huanglei/tools/gnuplot-install/bin/gnuplot < plot.sh

    echo "<table border=\"1\" cellpadding=\"3\" cellspacing=\"3\"> <td><img src=\"plot-1d-$NumTorsion.png\"></td> </table><br><br><br>" >> result-1D.html
    
    NumTorsion=`expr $NumTorsion + $NINC`
  else
    NumTorsion=`expr $NumTorsion - $NINC`
    echo "Number of torsion = $NumTorsion"
    exit 0
  fi
done

