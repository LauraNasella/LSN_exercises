CC = g++
CFLAGS = -Wall -O3 --std=c++11

main8.2.exe : main8.2.o random.o
	$(CC) random.o main8.2.o -o main8.2.exe
main8.2.o : main8.2.cpp
	$(CC) -c main8.2.cpp -o main8.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main8.2.exe seed.out
