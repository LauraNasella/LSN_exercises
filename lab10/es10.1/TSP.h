#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

#ifndef __TSP__
#define __TSP__

class City {
	private:
		double x, y;
		int pos; //crossover
	public:
		City(double X, double Y) : x(X), y(Y){} 
		City(const City & c)
        	: x(c.x),
          	  y(c.y){}
          	~City();
		double getx() const{ return x; }
    		double gety() const{ return y; }
		void set_pos(int & posiz);
		int get_pos() {return pos;}
		City& operator=(const City&);
};

//Individuo rappresenta un percorso tra le città 
class Individuo { 
  	private:
    		int N; //Numero di città
		City partenza;
		vector<City> percorso;
		vector<City> path_totale;
  	public:
		Individuo(int number, const City & inizio, vector<City>::iterator it1, vector<City>::iterator it2);
		Individuo(const Individuo & i)
		: N(i.N),
		  partenza(i.partenza),
		  percorso(i.percorso){}
		~Individuo();
		int get_N() const {return N;}
		City get_partenza() const {return partenza;}
		vector<City> get_percorso() const {return percorso;}
		void mutation_swap(Random & rnd);
		void mutation_shift(Random & rnd);
		void mutation_permutation(Random & rnd);	
		void mutation_inversion(Random & rnd);
		double distanza2(int i, int j) const;
		double L1() const; 
		bool check();
		void set_path_totale();
		vector<City> get_path_totale() const {return path_totale;}
		Individuo& operator=(const Individuo&);
};		

//Popolazione di individui
class Popolazione {
	private:
		int N_pop;
		vector<Individuo> pop;
	public:
		Popolazione(int number, vector<Individuo>::iterator it1, vector<Individuo>::iterator it2);	
		Popolazione(int number, Random & rnd, Individuo primo); //Creo la prima popolazione dal primo Individuo generato
		~Popolazione();
		vector<Individuo> get_pop() const {return pop;}
		void sort_L1(); 
		int selezione(Random & rnd);
		void crossover(Random & rnd, int k, int l); //tra due individui di pop
		void new_generazione(Random & rnd); //per modificare la vecchia generazione con le mutazioni e il crossover
		double mean_L1();
		void migrazione(double * new_x, double * new_y); //serve nell'esercitazione 10
};

#endif // __TSP__
