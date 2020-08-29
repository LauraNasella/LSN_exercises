#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <cfenv>
#include <climits>

using namespace std;
 
double error(double * medie, double * medie2, int i);

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	ofstream file_out;
	ifstream file_in;
	int M=10000000; //N_throws
	int N=100; //N blocchi
	int R= M/N; //tiri in ogni blocco
	double d = 2.0; //ho scelto io 2, distanza tra le linee
	double L = 1.0; //sempre scelto io, tale che d>L

	bool * vettore_throws= NULL;	
	double * pi = NULL;
	double * pi2 = NULL;
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

	vettore_throws = new bool[M];

/*Considero ad esempio un quadrato 20x20, con linea ogni 2 in verticale. Scelgo due numeri casuali x_c e y_c in [0,20) (ma posso cambiarli a piacere): li considero come il centro di una circonferenza di raggio r=L=1. Devo trovare un altro punto sulla circonferenza, non posso farlo invertendo l'equazione della circonferenza perchè avrei dei problemi nell'estrazione uniforme, quindi uso un metodo accept-reject. A questo punto devo capire se l'ago interseca una linea oppure no. Approssimo le y e y_c agli interi più vicini con la funzione round() e considero tutti i vari casi per capire se l'ago interseca o meno la linea.
*/
	double xmin = 0.0;
	double xmax = 20.0;

	for (int i=0; i<M; i++){
		
		//double x_c = xmin + (xmax - xmin) * rnd.Rannyu();
		double y_c = xmin + (xmax - xmin) * rnd.Rannyu();

		//Uso un metodo accept-reject
    	 	double r=0.;
     		double s=0.;
     		do{
       			r=rnd.Rannyu(-1,+1);
       			s=rnd.Rannyu(-1,+1);
     		}while(r*r+s*s>1 || (r==0 && s==0));//Genero due punti nella circonferenza di raggio r=L=1
	
     		//double cos = r/sqrt(r*r+s*s);
		double sin = s/sqrt(r*r+s*s);
		
		//double x = cos + x_c;
		double y = sin + y_c;

		int yy_c = round(y_c);
		int yy = round(y);
	
		if (yy_c == yy){ //cioè se vengono approssimati allo stesso intero
			if( (yy % 2)==0 ){  //i multipli di 2 sono dove c'è una linea
				if (y_c > y){
					if  ((y_c > yy) && (y < yy)){ //interseca
						vettore_throws[i] = true;
					}
					else { //non interseca
						vettore_throws[i] = false;
					}
				}
				else { //cioè se y_c < y	
					if  ((y > yy) && (y_c < yy)){ //interseca
						vettore_throws[i] = true;
					}
					else { //non interseca
						vettore_throws[i] = false;
					}
				}
			}
			else { //se sono uguali ma non multipli di 2, non interseca sicuramente, cioè l'ago è in mezzo a due linee
				vettore_throws[i]=false;	
			}	
		}
		else { // se non viene lo stesso intero ho 4 casi possibili. Sicuramente uno dei due è pari
			if((yy % 2)==0){
				if( ((y>yy) && (y_c>yy)) || ((y<yy) && (y_c<yy)) ) { //se sono entrambi maggiori o minori di quello pari non interseca
					vettore_throws[i] = false;
				}
				else{	
					vettore_throws[i] = true;
				}
			}
			else{ //vuole dire che yy_c è quello pari
				if( ((y>yy_c) && (y_c>yy_c)) || ((y<yy_c) && (y_c<yy_c)) ) { //se sono entrambi maggiori o minori di quello pari non interseca
					vettore_throws[i] = false;
				}
				else{	
					vettore_throws[i] = true;
				}
			}
		}		
	}

	rnd.SaveSeed();
 
	pi = new double[N];
	pi2 = new double[N];
	sum_prog = new double[N];
	sum2_prog = new double[N];
	err_prog = new double[N];
	x = new int[N];

	for (int i=0; i<N; i++){
    		double hit = 0;
   		for (int j=0; j<R; j++){
        		int k = j+i*R;
			if(vettore_throws[k]==true){ //vuol dire che l'ago aveva intersecato la linea
        			hit = hit + 1;
			}
			else{
			}
		}
		pi[i] = (2*L*R)/(hit*d);	//Il numero di tiri corrisponde a R=M/N, cioè quanti elementi ci sono nel blocco i-esimo che sto considerando degli N
    		pi2[i] = pow(pi[i],2);
	}

	
	for (int i=0; i<N; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + pi[j];
			sum2_prog[i] = sum2_prog[i]+ pi2[j];
	
		}
		sum_prog[i] = sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i] = sum2_prog[i]/(i+1); //Cumulative square average
   		err_prog[i] = error(sum_prog,sum2_prog,i); //Statistical uncertainty
	}
	
	file_out.open("pi_greco.out");	
	
	for(int i=0; i<N; i++){
		x[i]=i;
	}

	for(int i=0; i<N; i++){
		file_out << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
	}

	file_out.close();

	return 0;
}

double error(double * sum_prog, double * sum2_prog, int i){ 

	if (i==0){
	        return 0;
	}
	else{
        	return sqrt((sum2_prog[i] - pow(sum_prog[i],2))/(i-1));	
	 }
}

