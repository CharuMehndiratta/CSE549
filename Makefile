CC = g++
CFLAGS  = -g -Wall

all:
	g++ min_containment_hash.cpp  MurmurHash3.cpp -o min_containment_hash -I include/

clean:
	$(RM) min_containment_hash *.o *~
