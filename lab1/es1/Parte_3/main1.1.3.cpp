#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
double chi_quadro(int * n_i, int M, int n); //Funzione per calcolare il chi-quadro

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	int M=100; //divido [0,1] in 100 intervalli, utilizzo M=100 anche per indicare quante volte ripeto il calcolo del chi-quadro
	int n=10000; //prendo 10000 numeri pseudo-casuali ogni volta
	
	double * vettore_casuali= NULL;	
	double * chi = NULL;

	ofstream file_out;

	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
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
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	vettore_casuali = new double[n];
	chi = new double[M];

	for (int i=0; i<M; i++){ //ripeto M=100 volte con ogni volta n=10000 numeri casuali

		int * n_i= NULL; //Mi serve per contare quanti eventi ci sono in ognuno dei 100 bin in cui ho diviso [0,1]
		n_i = new int[M];
		
		for(int j=0; j<n; j++){
			vettore_casuali[j] = rnd.Rannyu();
			++n_i[int(vettore_casuali[j]*M)]; //così conto 
		}

		chi[i] = chi_quadro(n_i,M,n);
	}

	rnd.SaveSeed();
	
	file_out.open("chi_quadro.out");	
	
	for(int i=0; i<M; i++){
		file_out << i << " " << chi[i] << endl;
	}

	file_out.close();

	return 0;
}

double chi_quadro(int * n_i, int M, int n){ //Funzione per calcolare il chi-quadro

	double expected_value = n/M; 
	double chiquadro = 0;

	for (int i=0; i<M; i++){
		chiquadro = chiquadro + pow((n_i[i]-expected_value),2)/expected_value;
	}

return chiquadro;   
}

