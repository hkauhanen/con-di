#!/bin/sh

pdflatex eigenplot1
pdflatex eigenplot2
pdflatex eigenplot
pdflatex bifuplot
pdflatex sweepplot

cp bifuplot.pdf ../plots/Fig1.pdf
cp sweepplot.pdf ../plots/Fig2.pdf
cp eigenplot.pdf ../plots/Fig3.pdf
