#!/bin/bash

gcc -g -fopenmp mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -I. -lm -o mcml.out
