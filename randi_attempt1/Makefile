
CC = gcc
MPI = mpicc

all: main.c serialtester.c datatrim.c
	$(CC) serialtester.c Lab4_IO.c -o serialtester -lm
	$(CC) datatrim.c -o datatrim
	$(MPI) main.c Lab4_IO.c -o main -lm

clean:
	-rm main
	-rm serialtester
	-rm datatrim
