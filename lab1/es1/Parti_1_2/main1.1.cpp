#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
double error(double * medie, double * medie2, int i);

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	ofstream file_out;
	ifstream file_in;
	int M=10000;
	int N=100;
	int L= M/N;
	double * vettore_casuali= NULL;	
	double * medie = NULL;
	double * medie2 = NULL;
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

	vettore_casuali = new double[M];

	for (int i=0; i<M; i++){
		vettore_casuali[i] = rnd.Rannyu();
	}

	rnd.SaveSeed();

//1.1 Media  

	medie = new double[N];
	medie2 = new double[N];
	sum_prog = new double[N];
	sum2_prog = new double[N];
	err_prog = new double[N];
	x = new int[N];

	for (int i=0; i<N; i++){
    		double sum1 = 0;
   		for (int j=0; j<L; j++){
        		int k = j+i*L;
        		sum1 = sum1 + vettore_casuali[k];
			
		}
		medie[i] = sum1/L;     
    		medie2[i] = pow(medie[i],2);
	}

	for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + medie[j];
			sum2_prog[i] = sum2_prog[i]+ medie2[j];
	
		}
		sum_prog[i] = sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i] = sum2_prog[i]/(i+1); //Cumulative square average
   		err_prog[i] = error(sum_prog,sum2_prog,i); //Statistical uncertainty
	}
	
	file_out.open("medie.out");	
	
	for(int i=0; i<N; i++){
		x[i]=i;
	}

	for(int i=0; i<N; i++){
		file_out << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
	}

	file_out.close();

//1.2 Varianza

	ofstream file_out_var;

	double * var = NULL;
	double * var2 = NULL;
	double * sum_prog_var = NULL;
	double * sum2_prog_var = NULL;
	double * err_prog_var = NULL;

	var = new double[N];
	var2 = new double[N];
	sum_prog_var = new double[N];
	sum2_prog_var = new double[N];
	err_prog_var = new double[N];

	for (int i=0; i<N; i++){
    		double sum = 0;
   		for (int j=0; j<L; j++){
        		int k = j+i*L;
        		sum = sum + pow((vettore_casuali[k]-0.5),2);	
		}
		var[i] = sum/L; //stima della varianza in ogni blocco di 100 "misure"
    		var2[i] = pow(var[i],2);
	}

	for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog_var[i] = sum_prog_var[i] + var[j];
			sum2_prog_var[i] = sum2_prog_var[i]+ var2[j];
	
		}
		sum_prog_var[i] = sum_prog_var[i]/(i+1); //Cumulative average
    		sum2_prog_var[i] = sum2_prog_var[i]/(i+1); //Cumulative square average
   		err_prog_var[i] = error(sum_prog_var,sum2_prog_var,i); //Statistical uncertainty
	}
	
	file_out_var.open("varianze.out");	
	
	for(int i=0; i<N; i++){
		file_out_var << x[i] << " " << sum_prog_var[i] << " " << err_prog_var[i] << endl;
	}
	
	file_out_var.close();

	return 0;
}

double error(double * sum_prog, double * sum2_prog, int i){ //Function for statistical uncertainty estimation

	if (i==0){
	        return 0;
	}
	else{
        	return sqrt((sum2_prog[i] - pow(sum_prog[i],2))/(i-1));	
	 }
}

