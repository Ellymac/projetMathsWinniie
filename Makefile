CC=g++
FLAGS=-Wall -Wextra

all: src/mhd.cpp
	$(CC) src/mhd.cpp -o mhd $(FLAGS)
