CC=g++
FLAGS=-Wall -Wextra

resolve.o: src/resolve.cpp src/resolve.h
	$(CC) -c src/resolve.cpp $(FLAGS)

all: src/mhd.cpp resolve.o
	$(CC) src/mhd.cpp resolve.o -o mhd $(FLAGS)

clean:
	rm -rf mhd clmhd.msh *.o
