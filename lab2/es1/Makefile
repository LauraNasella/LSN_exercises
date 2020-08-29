CC = g++
CFLAGS = -Wall -O3 --std=c++11

main2.1.exe : main2.1.o random.o
	$(CC) random.o main2.1.o -o main2.1.exe
main2.1.o : main2.1.cpp
	$(CC) -c main2.1.cpp -o main2.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main2.1.exe seed.out
