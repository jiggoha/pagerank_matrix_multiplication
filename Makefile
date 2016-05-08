TIMINGDIR = /home/fas/hpcprog/ahs3/cpsc424/utils/timing
# CC = icc
# CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I$(TIMINGDIR)

matrix_multiplication: matrix_multiplication.o random_list.o $(TIMINGDIR)/timing.o
	mpicc -o $@ $(CFLAGS) $^

test: test.o random_list.o
	gcc test.c -o test -g

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f matrix_multiplication *.o 
