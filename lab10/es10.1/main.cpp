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

double Metropolis(double L_old, double L_new, double beta, int * n_acc, int * n_tot, Random rnd);

/*
Passo gli argomenti dalla riga di comando: 
-) 0 se voglio usare l'algoritmo genetico della 9
-) 1 se voglio usare il simulated annealing
*/

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


	int algoritmo = atoi(argv[1]);
		
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
	
if(algoritmo == 0){ //ALGORITMO GENETICO

	ofstream best_path_circ, best_half_circ, final_cities_circ;
	ofstream best_path_quad, best_half_quad, final_cities_quad;

	best_path_circ.open("L1_best_path_circ.out");
	best_half_circ.open("L1_best_half_circ.out");
	final_cities_circ.open("final_cities_circ.out");
	
	best_path_quad.open("L1_best_path_quad.out");
	best_half_quad.open("L1_best_half_quad.out");
	final_cities_quad.open("final_cities_quad.out");
	
	//primo individuo
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
	int n_generazioni = 600;	
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

	n_generazioni = 2000;
	
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
}

else{ //SIMULATED ANNEALING

	ofstream SA_best_path_circ, SA_final_cities_circ;
	ofstream SA_best_path_quad, SA_final_cities_quad;

	SA_best_path_circ.open("SA_L1_best_path_circ.out");
	SA_final_cities_circ.open("SA_final_cities_circ.out");
	
	SA_best_path_quad.open("SA_L1_best_path_quad.out");
	SA_final_cities_quad.open("SA_final_cities_quad.out");

	double beta_i = 0.02;
	double beta;
	double L_old = 0.;
	double L_old_q = 0.;
	double L_new = 0.;
	double L_new_q = 0.;
	int n_steps = 4000;
	int n_i = 100;
	int n_acc = 0;
    	int n_tot = 0;
    	int n_acc_q = 0;
    	int n_tot_q = 0;
    	
	//Primo individuo
	Individuo salesman(N, inizio, cities.begin(), cities.end());

	if (salesman.check()==false){ //false per me vuol dire che non passa nella stessa città
	}
	else{
		cout <<"Il percorso passa due volte nella stessa città! Non va bene!" << endl;
	}
	
	Individuo salesman_q(N, inizio_q, cities_q.begin(), cities_q.end());
	if (salesman_q.check()==false){ //false per me vuol dire che non passa nella stessa città
	}
	else{
		cout <<"Il percorso passa due volte nella stessa città! Non va bene!" << endl;
	}
	
	for(int i=0; i<n_steps; i++){
		beta = beta_i + i*0.025;
		for(int j=0; j<n_i; j++){
			L_old = salesman.L1();
			L_old_q = salesman_q.L1();
			Individuo copia(salesman);
			Individuo copia_q(salesman_q);
			
			int quale = rnd.Rannyu(0,4);
			if (quale==0) copia.mutation_swap(rnd);
      			else if (quale==1) copia.mutation_shift(rnd);
      			else if (quale==2) copia.mutation_permutation(rnd);
      			else copia.mutation_inversion(rnd); 
      			
      			quale = rnd.Rannyu(0,4);
			if (quale==0) copia_q.mutation_swap(rnd);
      			else if (quale==1) copia_q.mutation_shift(rnd);
      			else if (quale==2) copia_q.mutation_permutation(rnd);
      			else copia_q.mutation_inversion(rnd); 
      			
      			L_new = copia.L1();
      			L_new_q = copia_q.L1();
      			double L1 = Metropolis(L_old,L_new,beta,& n_acc, &n_tot, rnd);
      			double L1_q = Metropolis(L_old_q,L_new_q,beta,& n_acc_q, &n_tot_q, rnd);
      			if (L1==L_new){
				salesman = copia;
			}
			else{ //niente
			}
			if (L1_q==L_new_q){
				salesman_q = copia_q;
			}
			else{ //niente
			}
			SA_best_path_circ << L1 << endl;
			SA_best_path_quad << L1_q << endl;	
		}
	}
	
	SA_final_cities_circ << salesman.get_partenza().getx() << " " << salesman.get_partenza().gety() << endl;

	for (auto j : salesman.get_percorso()){
		SA_final_cities_circ << j.getx() << " " << j.gety() << endl;
	}
	
	SA_final_cities_quad << salesman_q.get_partenza().getx() << " " << salesman_q.get_partenza().gety() << endl;

	for (auto j : salesman_q.get_percorso()){
		SA_final_cities_quad << j.getx() << " " << j.gety() << endl;
	}
			
	SA_best_path_quad.close();
	SA_final_cities_quad.close();
	SA_best_path_circ.close();
	SA_final_cities_circ.close();
}
   	return 0;
}
double Metropolis(double L_old, double L_new, double beta, int * n_acc, int * n_tot, Random rnd){
	double rapporto = exp(-beta*(L_new-L_old));
	
	double alpha = 0.;
	
	if(rapporto < 1){
		alpha = rapporto;
	}
	else{
		alpha = 1;
	}
	
	double r = rnd.Rannyu(); //genero un numero casuale uniformemente tra 0 e 1
	if(r<=alpha){
		*n_acc = *n_acc +1;
		*n_tot = *n_tot +1;
		return L_new;
	}	
	else{
		*n_tot = * n_tot +1;
		return L_old;
	}		
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
