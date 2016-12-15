CC=gcc
CFLAGS=-O3
LFLAGS=-lm -lgsl -lcblas

run: f1
	./f1

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

f1: main.o f1.o
	$(CC) $(LFLAGS) -o $@ f1.o main.o

clean:
	rm f1 f2 *.o
