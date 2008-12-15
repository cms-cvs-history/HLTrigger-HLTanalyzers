#!/bin/bash

#echo $1
latex $1.tex
latex $1.tex
dvipdf $1.dvi
rm $1.dvi
rm $1.aux
rm $1.log
#rm $1*~
acroread $1.pdf &
