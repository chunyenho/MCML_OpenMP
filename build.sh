#!/bin/bash

gcc -fopenmp mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -I. -lm -o mcml.out
