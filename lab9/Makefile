CC = g++
CCC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o TSP.o
	$(CCC) random.o main.o TSP.o -o main.exe
main.o : main.cpp TSP.h random.h
	$(CCC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
TSP.o : TSP.cpp TSP.h random.h
	$(CC) -c TSP.cpp -o TSP.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
