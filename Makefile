CC=gcc
CFLAGS=-Wall -W -pedantic
CLIB=-lm

moco: main.o motif.o stat.o
	$(CC) *.o $(CLIB) -o moco
main.o: main.c
	$(CC) $(CFLAGS) main.c -c
motif.o: motif.c motif.h
	$(CC) $(CFLAGS) motif.c -c
stat.o: stat.c stat.h
	$(CC) $(CFLAGS) stat.c -c
clean:
	rm -f *.o moco
