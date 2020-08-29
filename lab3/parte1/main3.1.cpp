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
	ofstream file_out1;
	//ifstream file_in;

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


	int M=100000; 
    	int N=100; //numero di blocchi
    	int L=int(M/N); //numero di stime in ogni blocco
    
	double S0=100.;
    	double T=1.;
    	double K=100.;
    	double r=0.1;
    	double sigma=0.25;
    	
	double * call = NULL;
	double * call2 = NULL;
	double * sum_prog_call = NULL;
	double * sum2_prog_call = NULL;
	double * err_prog_call = NULL;
	double * put = NULL;
	double * put2 = NULL;
	double * sum_prog_put = NULL;
	double * sum2_prog_put = NULL;
	double * err_prog_put = NULL;
	int * x = NULL;
	
//1)Campionamento diretto
      	
      	call = new double[N];
      	call2 = new double[N];
      	put = new double[N];
      	put2 = new double[N];
	sum_prog_call = new double[N];
	sum2_prog_call = new double[N];
	err_prog_call = new double[N];
	sum_prog_put = new double[N];
	sum2_prog_put = new double[N];
	err_prog_put = new double[N];
	x = new int[N];
     
    	for(int i=0; i<N; i++){
        	double c = 0.;
        	double p = 0.;
        	for(int j=0; j<L; j++){
            		double w = rnd.Gauss(0,T);
            		double ST = S0 * exp((r - 0.5 * pow(sigma,2)) * T + sigma * w);
            		c = c + exp(-r*T) * max(0.,ST-K);
            		p = p + exp(-r*T) * max(0.,K-ST);
        	}
        	call[i] = c/L;
        	call2[i] = pow(call[i],2);
        	put[i] = p/L;
        	put2[i] = pow(put[i],2);	
        }
        
        for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog_call[i] = sum_prog_call[i] + call[j];
			sum2_prog_call[i] = sum2_prog_call[i]+ call2[j];
			sum_prog_put[i] = sum_prog_put[i] + put[j];
			sum2_prog_put[i] = sum2_prog_put[i]+ put2[j];
		}
		sum_prog_call[i] = sum_prog_call[i]/(i+1); //Cumulative average
    		sum2_prog_call[i] = sum2_prog_call[i]/(i+1); //Cumulative square average
   		err_prog_call[i] = error(sum_prog_call,sum2_prog_call,i); //Statistical uncertainty
   		
   		sum_prog_put[i] = sum_prog_put[i]/(i+1); //Cumulative average
    		sum2_prog_put[i] = sum2_prog_put[i]/(i+1); //Cumulative square average
   		err_prog_put[i] = error(sum_prog_put,sum2_prog_put,i); //Statistical uncertainty
	}
	
	file_out.open("Call_diretto.out");
   	file_out1.open("Put_diretto.out");	
	
	for(int i=0; i<N; i++){
		x[i]=i;
	}

	for(int i=0; i<N; i++){
		file_out << x[i] << " " << sum_prog_call[i] << " " << err_prog_call[i] << endl;
		file_out1 << x[i] << " " << sum_prog_put[i] << " " << err_prog_put[i] << endl;
	}

	file_out.close();
	file_out1.close();

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
