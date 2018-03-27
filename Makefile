
CC1 = gcc
CC2 = mpicc

all: main.c serialtester.c datatrim.c
	$(CC1) serialtester.c Lab4_IO.c -o serialtester -lm
	$(CC1) datatrim.c -o datatrim
	$(CC2) -g -Wall main.c -o main -lm

clean:
	-rm main
	-rm serialtester
	-rm datatrim
