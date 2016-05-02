#!/bin/bash

icc -O3  -vec-report2 mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -fopenmp -I. -lm -o mcml.out

