CC = gcc
FLAGS = -fopenmp -O3 -vec-report2 
MIC_FLAGS = -fopenmp -O3 -vec-report2 -mmic
INCLUDE = .



all:
	$(CC) $(FLAGS) mcmlmain.c mcmlio.c mcmlnr.c mcmlgo.c -I$(INCLUDE) -o mcml_xeon.out
	$(CC) $(MIC_FLAGS) mcmlmain.c mcmlio.c mcmlnr.c mcmlgo.c -I$(INCLUDE) -o mcml_xphi.out
 
