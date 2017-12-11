CC = g++
CFLAGS  = -g -Wall

all:
	g++ query.cpp MurmurHash3.cpp -o query -I include/
	g++ min_containment_hash.cpp MurmurHash3.cpp -o min_containment_hash -I include/

clean:
	$(RM) query *.o *~
