CC = g++
CFLAGS  = -g -Wall

all:
	g++ min_containment_hash.cpp MurmurHash3.cpp utils.cpp -o min_containment_hash -I include/
	g++ query.cpp MurmurHash3.cpp utils.cpp -o query -I include/

clean:
	$(RM) min_containment_hash *.o *~
	$(RM) query *.o *~
