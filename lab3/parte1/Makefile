CC = g++
CFLAGS = -Wall -O3 --std=c++11

main3.1.exe : main3.1.o random.o
	$(CC) random.o main3.1.o -o main3.1.exe
main3.1.o : main3.1.cpp
	$(CC) -c main3.1.cpp -o main3.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main3.1.exe seed.out
