CC		= gcc 
WARNS	= -Wextra -Wno-unused-parameter -Wall -Wformat=2 -Wuninitialized \
		  -Wfloat-equal -Wshadow -Wpointer-arith \
		  -Wstrict-prototypes -Wmissing-prototypes -Wmaybe-uninitialized

REAL	= double
		  
LFLAGS	= -std=c99 -pedantic -O3 -fPIC $(WARNS) -Werror -DREAL=$(REAL)
TFLAGS	= -std=c99 -pedantic -O3 $(WARNS) -Werror -DREAL=$(REAL)

LIBS	= -lm

OBJ		= loci.o 1D.o 2D.o 3D.o 4D.o 
TOBJ	= runtests.o 1Dtest.o 2Dtest.o 3Dtest.o 4Dtest.o

BIN		= bin
SRC		= src
INCL	= include
TEST	= test

PYTEST	= python/test

vpath %.c  $(SRC):$(TEST)
vpath %.h  $(INCL):$(TEST):$(SRC)

shlib :	$(OBJ)
		$(CC) $(LFLAGS) -shared $(OBJ) -o lib/libloci.so

test :	tfiles
		$(BIN)/runtests

pytest : shlib
		python $(PYTEST)/runtests.py

tfiles : $(OBJ) $(TOBJ)
		$(CC) $(TFLAGS) $(LIBS) $(OBJ) $(TOBJ) -o $(BIN)/runtests

runtests.o :	runtests.c runtests.h 
		$(CC) $(CFLAGS) $(LIBS) -c $<

loci.o :	loci.c loci.h 
		$(CC) $(LFLAGS) $(LIBS) -c $<

1D.o :	1D.c 1D.h  loci.h helpers.h
		$(CC) $(LFLAGS) $(LIBS) -c $< 

2D.o :	2D.c 2D.h  loci.h helpers.h
		$(CC) $(LFLAGS) $(LIBS) -c $< 

3D.o :	3D.c 3D.h  loci.h helpers.h
		$(CC) $(LFLAGS) $(LIBS) -c $< 

4D.o :	4D.c 4D.h  loci.h helpers.h
		$(CC) $(LFLAGS) $(LIBS) -c $< 

1Dtest.o : 1Dtest.c test/runtests.h test/testhelpers.h
		$(CC) $(TFLAGS) $(LIBS) -c $< 

2Dtest.o : 2Dtest.c test/runtests.h test/testhelpers.h
		$(CC) $(TFLAGS) $(LIBS) -c $< 

3Dtest.o : 3Dtest.c test/runtests.h test/testhelpers.h
		$(CC) $(TFLAGS) $(LIBS) -c $< 

4Dtest.o : 4Dtest.c test/runtests.h test/testhelpers.h
		$(CC) $(TFLAGS) $(LIBS) -c $< 

.PHONY : clean
clean :	
		rm -f *.o lib/* $(BIN)/*
