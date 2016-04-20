#!/bin/bash

gcc -fopenmp -std=c99 mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -I. -lm -o mcml.out
 
