#!/bin/bash

cat "script_plot_1.txt" > plot.plt

./Debug/Genetisch >> plot.plt

cat "script_plot_2.txt" >> plot.plt

gnuplot "plot.plt"
