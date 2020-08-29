CC = g++
CFLAGS = -Wall -O3 --std=c++11

main1.1.3.exe : main1.1.3.o random.o
	$(CC) random.o main1.1.3.o -o main1.1.3.exe
main1.1.3.o : main1.1.3.cpp
	$(CC) -c main1.1.3.cpp -o main1.1.3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main1.1.3.exe seed.out
