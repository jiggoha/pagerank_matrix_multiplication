TIMINGDIR = /home/fas/hpcprog/ahs3/cpsc424/utils/timing
CC = gcc
CFLAGS = -g -O0 -std=c99 -I$(TIMINGDIR)
# -xHost -fno-alias
# $(TIMINGDIR)/timing.o   <--- add to matrix_multiplication line when actually working on omega

matrix_multiplication: matrix_multiplication.o random_list.o 
	$(CC) -o $@ $(CFLAGS) $^

test: test.o random_list.o
	gcc test.c -o test -g

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f matrix_multiplication *.o 
