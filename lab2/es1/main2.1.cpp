#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
double error(double * medie, double * medie2, int i);

double g_x(double x);

double d_funzione(double x); 

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	ofstream file_out;
	int M=10000; //numero di lanci
	int N=100; //numero di blocchi
	int L= M/N; //numero di lanci in ogni blocco
	double * vettore_casuali= NULL;	
	double * I = NULL;
	double * I2 = NULL;
	double * sum_prog = NULL;
	double * sum2_prog = NULL;
	double * err_prog = NULL;
	int * x = NULL;

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

//1)Campionamento uniforme  

	vettore_casuali = new double[M];

	for (int i=0; i<M; i++){
		vettore_casuali[i] = rnd.Rannyu();
	}

	rnd.SaveSeed();

	I = new double[N];
	I2 = new double[N];
	sum_prog = new double[N];
	sum2_prog = new double[N];
	err_prog = new double[N];
	x = new int[N];

	for (int i=0; i<N; i++){
    		double sum1 = 0;
		double g = 0;
   		for (int j=0; j<L; j++){
        		int k = j+i*L;
			g = g_x(vettore_casuali[k]);
			sum1 = sum1 + g;
		}
		I[i] = sum1/L;     
    		I2[i] = pow(I[i],2);
	}

	for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + I[j];
			sum2_prog[i] = sum2_prog[i]+ I2[j];
	
		}
		sum_prog[i] = sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i] = sum2_prog[i]/(i+1); //Cumulative square average
   		err_prog[i] = error(sum_prog,sum2_prog,i); //Statistical uncertainty
	}
	
	file_out.open("I_uniforme.out");	
	
	for(int i=0; i<N; i++){
		x[i]=i;
	}

	for(int i=0; i<N; i++){
		file_out << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
	}

	file_out.close();

// Importance sampling

	for (int i=0; i<M; i++){
		vettore_casuali[i] = rnd.d_x(); //lo riempio ora secondo la nuova distribuzione d
	}

	rnd.SaveSeed();

	for (int i=0; i<N; i++){
    		double sum1 = 0;
		double g = 0;
		double d = 0;
   		for (int j=0; j<L; j++){
        		int k = j+i*L;
			g = g_x(vettore_casuali[k]);
			d = d_funzione(vettore_casuali[k]);
			sum1 = sum1 + g/d;
		}
		I[i] = sum1/L;     
    		I2[i] = pow(I[i],2);
	}

	for (int i=0; i<N; i++){
		sum_prog[i]=0;
		sum2_prog[i]=0;
	}

	for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + I[j];
			sum2_prog[i] = sum2_prog[i]+ I2[j];
	
		}
		sum_prog[i] = sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i] = sum2_prog[i]/(i+1); //Cumulative square average
   		err_prog[i] = error(sum_prog,sum2_prog,i); //Statistical uncertainty
	}
	
	file_out.open("I_importance.out");	
	
	for(int i=0; i<N; i++){
		x[i]=i;
	}

	for(int i=0; i<N; i++){
		file_out << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
	}

	file_out.close();

	return 0;
}

double error(double * medie, double * medie2, int i){ 

	if (i==0){
	        return 0;
	}
	else{
        	return sqrt((medie2[i] - pow(medie[i],2))/(i-1));	
	 }
}

double g_x(double x){
	return ((M_PI)/2)*cos((M_PI*x)/2);
}

double d_funzione(double x){
	return -2*x+2;
}

