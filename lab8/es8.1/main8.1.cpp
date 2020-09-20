#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
double error(double * medie, double * medie2, int i);

double psi_trial (double x, double mu, double sigma);

void Metropolis(double * pos0, double * pos_new, int * n_acc, int * n_tot, Random rnd, double mu, double sigma);

double V(double x){ return pow(x,4) -2.5*pow(x,2);}

double H(double x, double mu, double sigma);
//Passo gli argomenti dalla riga di comando: dopo il nome del programma passo mu e sigma

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

	if( argc != 3 ){
		cerr << "ATTENZIONE: il programma ha bisogno di 3 argomenti, il primo è il nome del programma e poi mu e sigma" << endl;
    		exit(1);
    	}
    	
    	double mu = atof(argv[1]);
    	double sigma = atof(argv[2]);

	int M=100000;
	int M_eq = 1000; //aggiungo 1000 steps che considero come equilibrazione
    	int N=100; //numero di blocchi
    	int L=int(M/N); //numero di stime in ogni blocco
    	
    	file_out.open("ave_H.out");
    	file_out1.open("campionamento.out");

    	//Scelgo il punto di partenza, devo sceglierlo in modo ragionevole
    	double pos0 = 0;
    	
    	//Devo calcolare l'ampiezza dello step che mi dia un'accettazione di ~50%.	
    	int n_acc = 0;
    	int n_tot = 0;
    	
    	double step = 4.5*sigma;
    	
    	double * H_vec = NULL;
    	double * H_mean = NULL;
	double * H_mean2 = NULL;
	double * sum_prog = NULL;
	double * sum2_prog = NULL;
	double * err_prog = NULL;
	int * x = NULL;
	
      	H_vec = new double[M];
      	H_mean = new double[N];
      	H_mean2 = new double[N];
	sum_prog = new double[N];
	sum2_prog = new double[N];
	err_prog = new double[N];
	x = new int[N];
	
    	for(int i=0; i<(M+M_eq); i++){
    		double pos_new = 0;
    		// trasizione uniforme
    		pos_new = pos0 + step * rnd.Rannyu(-1,+1);
    	
   		//Ora devo decidere se accettare la nuova posizione con Metropolis
   		Metropolis(& pos0, & pos_new, & n_acc, & n_tot, rnd, mu, sigma);
   
   		if (i%1000==0){
   			cout << "Accettazione: " <<  n_acc/double(n_tot) << endl;
   		}
   		
   		if(i>=M_eq){ //dopo equilibrazione
   			file_out1 << pos0 << endl;
   			file_out1 << -pos0 << endl;
   			H_vec[i-M_eq] = H(pos0,mu,sigma);
   			//cout << pos0 << " " << H_vec[i-M_eq] << endl;
   		}
    	}
    	  
    	for(int i=0; i<N; i++){
        	double sum = 0.;
        	for(int j=0; j<L; j++){
        		int k = j + i*L;
            		sum = sum + H_vec[k];
            		//cout << H_vec[k] << endl;
        	}
        	H_mean[i] = sum/L;
        	H_mean2[i] = pow(H_mean[i],2);	
        }
        
        for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + H_mean[j];
			sum2_prog[i] = sum2_prog[i]+ H_mean2[j];
		}
		sum_prog[i] = sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i] = sum2_prog[i]/(i+1); //Cumulative square average
   		err_prog[i] = error(sum_prog,sum2_prog,i); //Statistical uncertainty
	}
	
	for(int i=0; i<N; i++){
		x[i]=i;
	}

	for(int i=0; i<N; i++){
		file_out << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
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

double psi_trial (double x, double mu, double sigma) {
	return exp(-0.5*pow((x-mu)/sigma,2)) + exp(-0.5*pow((x+mu)/sigma,2));
}

void Metropolis(double * pos0, double *pos_new, int * n_acc, int * n_tot, Random rnd, double mu, double sigma){
	
	double rapporto = 0.;
	rapporto = pow(psi_trial(*pos_new,mu,sigma),2)/pow(psi_trial(*pos0,mu,sigma),2);

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
		*pos0 = *pos_new;
	}	
	else{
		*n_tot = * n_tot +1;
	}		
}

double H(double x, double mu, double sigma){

	double g_meno = exp(-0.5*pow((x-mu)/sigma,2));		
	double g_piu = exp(-0.5*pow((x+mu)/sigma,2));		

	double H_psi = 0.5*pow(sigma, -2) * ( (1 - pow((x-mu)/sigma , 2))*g_meno + (1 - pow((x+mu)/sigma , 2))*g_piu ); 

	return H_psi/psi_trial(x, mu, sigma) + V(x);
}

