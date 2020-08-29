CC = g++
CFLAGS = -Wall -O3 --std=c++11

main2.2.exe : main2.2.o random.o
	$(CC) random.o main2.2.o -o main2.2.exe
main2.2.o : main2.2.cpp
	$(CC) -c main2.2.cpp -o main2.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main2.2.exe seed.out
