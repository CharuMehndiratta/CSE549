CC = g++
CFLAGS  = -g -Wall

all:
	g++ sequence_similarity.cpp  MurmurHash3.cpp -o sequence_similarity -I include/

clean:
	$(RM) sequence_similarity *.o *~
