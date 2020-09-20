CC = g++
CFLAGS = -Wall -O3 --std=c++11

main8.1.exe : main8.1.o random.o
	$(CC) random.o main8.1.o -o main8.1.exe
main8.1.o : main8.1.cpp
	$(CC) -c main8.1.cpp -o main8.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main8.1.exe seed.out
