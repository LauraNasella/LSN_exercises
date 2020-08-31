#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
double error(double * medie, double * medie2, int i);

double psi_100(double * pos);

double psi_210(double * pos);

void Metropolis(double * pos0, double * pos_new, int * n_acc, int * n_tot, int stato, Random rnd);
/*
Passo gli argomenti dalla riga di comando: 
-) 0 o 1 in base a se voglio campionare il modulo quadro di psi_100 (ground state) o di psi_210 (primo stato eccitato);
-) 0 o 1 in base a se voglio usare una probabilità di transizione uniforme o gaussiana.
Quindi i comandi possibili sono:
-) ./main5.1.exe 0 0 per psi_100 con prob. di transizione uniforme
-) ./main5.1.exe 0 1 per psi_100 con prob. di transizione gaussiana
-) ./main5.1.exe 1 0 per psi_210 con prob. di transizione uniforme
-) ./main5.1.exe 1 1 per psi_210 con prob. di transizione gaussiana
*/
int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	ofstream file_out;
	ofstream file_out1;
	ofstream file_out2;
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
		cerr << "ATTENZIONE: il programma ha bisogno di 3 argomenti, il primo è il nome del programma e poi 0 o 1 in base a quale funzione d'onda si vuole campionare e con che prob. di transizione." << endl;
    		exit(1);
    	}
    	
    	int stato = atoi(argv[1]);
    	int trans = atoi(argv[2]);
    	
    	//cout << stato << endl;
    	if ((stato != 0) && (stato != 1)){
    		cerr << "ATTENZIONE: il programma può solo campionare il modulo quadro di psi_100 (che corrisponde a 0) o di psi_210 (che corrisponde a 1). Quindi reinserisci gli input" << endl;
    		exit(1);
    	}
    	if ((trans != 0) && (trans != 1)){
    		cerr << "ATTENZIONE: il programma può solo usare o la prob. di transizione uniforme (che corrisponde a 0) o quella gaussiana (che corrisponde a 1). Quindi reinserisci gli input" << endl;
    		exit(1);
    	}

	int M=100000;
	int M_eq = 1000; //aggiungo 1000 steps che considero come equilibrazione
    	int N=100; //numero di blocchi
    	int L=int(M/N); //numero di stime in ogni blocco
    	
    	if (stato == 0 && trans == 0){
    		file_out.open("100_uniforme.out");
    		file_out1.open("100_unif_eq.out");
    		file_out2.open("100_unif_3D");
    	}
    	else if (stato == 0 && trans == 1){
    		file_out.open("100_gaussiana.out");
    		file_out1.open("100_gauss_eq.out");
    		file_out2.open("100_gauss_3D");
    	}
    	else if (stato == 1 && trans == 0){
    		file_out.open("210_uniforme.out");
    		file_out1.open("210_unif_eq.out");
    		file_out2.open("210_unif_3D");
    	}
    	else if (stato == 1 && trans == 1){
    		file_out.open("210_gaussiana.out");
    		file_out1.open("210_gauss_eq.out");
    		file_out2.open("210_gauss_3D");
    	}
    	
    	//Scelgo il punto di partenza nello spazio 3D, devo sceglierlo in modo ragionevole
    	double * pos0 = new double[3];
    	if(stato == 0){
    		pos0[0] = 0.;
    		pos0[1] = 0.;
    		pos0[2] = 0.;
    	}
    	else{ //cioè stato eccitato 210
    		pos0[0] = 0.;
    		pos0[1] = 0.;
    		pos0[2] = 5.;
    	}
    	
    	//Devo calcolare l'ampiezza dello step che mi dia un'accettazione di ~50% in tutti e 4 i casi.	
    	int n_acc = 0;
    	int n_tot = 0;
    	
    	double step = 0.;
    	if (stato == 0 && trans == 0){
		step = 1.2;
    	}
    	else if (stato == 0 && trans == 1){
	    	step = 0.75;
    	}
    	else if (stato == 1 && trans == 0){
		step = 2.9;
    	}
    	else if (stato == 1 && trans == 1){
		step = 1.8;
    	}
    	
    	double * r = NULL;
    	double * r_eq = NULL;
    	double * r_mean = NULL;
	double * r_mean2 = NULL;
	double * sum_prog = NULL;
	double * sum2_prog = NULL;
	double * err_prog = NULL;
	int * x = NULL;
	
      	r = new double[M];
      	r_eq = new double[M_eq];
      	r_mean = new double[N];
      	r_mean2 = new double[N];
	sum_prog = new double[N];
	sum2_prog = new double[N];
	err_prog = new double[N];
	x = new int[N];
	
    	for(int i=0; i<(M+M_eq); i++){
    		double * pos_new = new double[3];
    		if(trans == 0){ //cioè uniforme
    			for(int j=0; j<3; j++){
    				pos_new[j] = pos0[j] + step * rnd.Rannyu(-1,+1);
    			}
    		}
    		else{ //cioè gaussiana
    			for(int j=0; j<3; j++){
    				pos_new[j] = pos0[j] + rnd.Gauss(0,step);
    			}
    		}
    	
   		//Ora devo decidere se accettare la nuova posizione con Metropolis
   		Metropolis(pos0, pos_new, & n_acc, & n_tot, stato, rnd);
   
   		if (i%1000==0){
   			cout << "Accettazione: " <<  n_acc/double(n_tot) << endl;
   		}
   		
   		//Salvo i primi 1000 steps per l'equilibrazione
   		if(i<1000){
   			r_eq[i] = sqrt(pos0[0]*pos0[0] + pos0[1]*pos0[1] + pos0[2]*pos0[2]);
   			file_out1 << r_eq[i] << endl;
   		}
   		else{
   			file_out2 << pos0[0] << " " << pos0[1] << " " << pos0[2] << endl;
   			r[i-1000] = sqrt(pos0[0]*pos0[0] + pos0[1]*pos0[1] + pos0[2]*pos0[2]);
   		}
    	}
    	  
    	for(int i=0; i<N; i++){
        	double sum = 0.;
        	for(int j=0; j<L; j++){
        		int k = j + i*L;
            		sum = sum + r[k];
        	}
        	r_mean[i] = sum/L;
        	r_mean2[i] = pow(r_mean[i],2);	
        }
        
        for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + r_mean[j];
			sum2_prog[i] = sum2_prog[i]+ r_mean2[j];
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
	file_out2.close();

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

double psi_100(double * pos){ //considero a_0 = 1, considero già il modulo quadro della fne d'onda. pos[0] = x, pos[1] = y, pos[2] = z.
	double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
 	return (1/M_PI) * exp(-2*r);
}

double psi_210(double * pos){ //considero a_0 = 1, considero già il modulo quadro della fne d'onda.
	double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
	double cos_theta = pos[2]/r;
	return ((r*r)*exp(-r)*cos_theta*cos_theta)/(32*M_PI);
}

void Metropolis(double * pos0, double * pos_new, int * n_acc, int * n_tot, int stato, Random rnd){
	double rapporto = 0.;
	
	if (stato==0){
		rapporto = psi_100(pos_new)/psi_100(pos0);
	}
	else{
		rapporto = psi_210(pos_new)/psi_210(pos0);
	} 
	
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
		for(int k=0; k<3; k++){
			pos0[k] = pos_new[k];
		}
	}	
	else{
		*n_tot = * n_tot +1;
	}		
}

