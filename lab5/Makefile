CC = g++
CFLAGS = -Wall -O3 --std=c++11

main5.1.exe : main5.1.o random.o
	$(CC) random.o main5.1.o -o main5.1.exe
main5.1.o : main5.1.cpp
	$(CC) -c main5.1.cpp -o main5.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main5.1.exe seed.out
