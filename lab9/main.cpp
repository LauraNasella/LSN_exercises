/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

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
#include "TSP.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes, questo" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in, questo" << endl;

	ofstream best_path_circ, best_half_circ, final_cities_circ;
	ofstream best_path_quad, best_half_quad, final_cities_quad;

	best_path_circ.open("L1_best_path_circ.out");
	best_half_circ.open("L1_best_half_circ.out");
	final_cities_circ.open("final_cities_circ.out");
	
	best_path_quad.open("L1_best_path_quad.out");
	best_half_quad.open("L1_best_half_quad.out");
	final_cities_quad.open("final_cities_quad.out");

//CIRCONFERENZA
//Creo casualmente le 32 città sulla circonferenza.
	int N = 32;
	int Npop = 100;
	double theta = rnd.Rannyu(0,2*M_PI);
	double s = sin(theta);
	double r = cos(theta);

//Prima città, la fisso
	City inizio(r,s);
	
//Le altre 31 città
	std::vector<City> cities;
	for (int i=0; i<N-1; i++){
		theta = rnd.Rannyu(0,2*M_PI);
		s = sin(theta);
		r = cos(theta);
		City c(r,s);
		cities.push_back(c);
	}
	
//Creo il primo individuo
	Individuo salesman(N, inizio, cities.begin(), cities.end());
	if (salesman.check()==false){ //false per me vuol dire che non passa nella stessa città
	}
	else{
		cout <<"Il percorso passa due volte nella stessa città! Non va bene!" << endl;
	}

//Creo la prima popolazione e la ordino in base a L1
	Popolazione prima(Npop, rnd, salesman); 
	prima.sort_L1();

//Generazioni
	int n_generazioni = 400;	
	for(int g=0; g<n_generazioni; g++){
		prima.new_generazione(rnd);
		prima.sort_L1();
		best_path_circ << g << " " << prima.get_pop()[0].L1() << endl;
		best_half_circ << g << " " << prima.mean_L1() << endl;
	}
	
	final_cities_circ << prima.get_pop()[0].get_partenza().getx() << " " << prima.get_pop()[0].get_partenza().gety() << endl;

	for (auto j : prima.get_pop()[0].get_percorso()){
		final_cities_circ << j.getx() << " " << j.gety() << endl;
	}
	/*
	for (auto i : prima.get_pop()){ 			
		cout << i.L1() << endl;
		}
	*/
	best_path_circ.close();
	best_half_circ.close();
	final_cities_circ.close();

	rnd.SaveSeed();
//QUADRATO

//Prima città, la fisso
	double x=rnd.Rannyu(-1,+1);
     	double y=rnd.Rannyu(-1,+1);
     	
	City inizio_q(x,y);

//Le altre 31 città
	std::vector<City> cities_q;
	for (int i=0; i<N-1; i++){
		x=rnd.Rannyu(-1,+1);
       	y=rnd.Rannyu(-1,+1);
		City c_q(x,y);
		cities_q.push_back(c_q);
	}
	
//Primo individuo
	Individuo salesman_q(N, inizio_q, cities_q.begin(), cities_q.end());
	if (salesman_q.check()==false){ //false per me vuol dire che non passa nella stessa città
	}
	else{
		cout <<"Il percorso passa due volte nella stessa città! Non va bene!" << endl;
	}
	
//Prima popolazione
	Popolazione prima_q(Npop, rnd, salesman_q);
	prima_q.sort_L1();

	n_generazioni = 400;
	
	for(int g=0; g<n_generazioni; g++){
		prima_q.new_generazione(rnd);
		prima_q.sort_L1();
		best_path_quad << g << " " << prima_q.get_pop()[0].L1() << endl;
		best_half_quad << g << " " << prima_q.mean_L1() << endl;
	}
	
	final_cities_quad << prima_q.get_pop()[0].get_partenza().getx() << " " << prima_q.get_pop()[0].get_partenza().gety() << endl;

	for (auto j : prima_q.get_pop()[0].get_percorso()){
		final_cities_quad << j.getx() << " " << j.gety() << endl;
	}
			
	best_path_quad.close();
	best_half_quad.close();
	final_cities_quad.close();

   	return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
