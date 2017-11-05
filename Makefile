CC = g++
CFLAGS  = -g -Wall

all:
	g++ main.cpp  MurmurHash3.cpp -o main -I include/

clean:
	$(RM) main *.o *~
