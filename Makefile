TIMINGDIR = /home/fas/hpcprog/ahs3/cpsc424/utils/timing
CC = icc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I$(TIMINGDIR)

matrix_multiplication: matrix_multiplication.o $(TIMINGDIR)/timing.o
	mpicc -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f matrix_multiplication *.o 
