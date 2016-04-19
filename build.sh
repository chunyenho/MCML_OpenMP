#!/bin/bash

gcc mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -fopenmp -I. -lm -o mcml.out
