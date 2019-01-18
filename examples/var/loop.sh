#!/bin/bash


#V=17Q1.0
#V=17Q2.0
#source /opt/conda/bin/activate analysis-$V 

V=2018-1.0
source /opt/conda/bin/activate analysis-$V

while [ 1 ]
do

    ipython runXS.py
    sleep 10

done



