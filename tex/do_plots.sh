#!/bin/sh

pdflatex eigenplot1
pdflatex eigenplot2
pdflatex eigenplot
pdflatex bifuplot
pdflatex sweepplot
pdflatex evalplot

cp bifuplot.pdf ../plots/Fig1.pdf
cp sweepplot.pdf ../plots/Fig2.pdf
cp evalplot.pdf ../plots/Fig3.pdf
cp eigenplot.pdf ../plots/Fig4.pdf
