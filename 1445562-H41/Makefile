
CC = gcc
MPI = mpicc

all: main_attempt1.c main_attempt2.c serialtester.c datatrim.c
	$(CC) serialtester.c Lab4_IO.c -o serialtester -lm
	$(CC) datatrim.c -o datatrim
	$(MPI) main_attempt1.c Lab4_IO.c -o main_attempt1 -lm
	$(MPI) main_attempt2.c Lab4_IO.c -o main -lm

clean:
	-rm main_attempt1
	-rm main
	-rm serialtester
	-rm datatrim
