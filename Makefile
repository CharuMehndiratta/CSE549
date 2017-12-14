CC = g++
CFLAGS  = -g -Wall

all:
	g++ -std=c++11 query.cpp MurmurHash3.cpp -o query -I include/
	g++ -std=c++11 min_containment_hash.cpp MurmurHash3.cpp -o build_indices -I include/

clean:
	$(RM) query build_indices *.o *~
