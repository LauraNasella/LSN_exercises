#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>

using namespace std;
 
double error(double medie[][100], double medie2[][100], int j, int i);

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	ofstream file_out;
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

//2) RW nel continuo
	int a = 1;
	int N=10000; //numero di RWs
	int n_steps = 100; //steps in ogni RW
	int n_blocchi = 100;
	int L = int(N/n_blocchi); //numero di RWs in ogni blocco

	double * vettore_r_N2 = NULL; //è il vettore di (r_n)^2 a ogni passo i [1,100]	
	double * vettore_radici = NULL; //è il vettore della radice delle medie finali per ogni passo i [1,100]	
	double * vettore_errori = NULL; //è il vettore delgli errori finali per ogni passo i [1,100]	
	double * x = NULL;
	double * y = NULL;
	double * z = NULL;
	int * n = NULL;
	double sum_prog[n_blocchi][100];
	double sum2_prog[n_blocchi][100];
	double err_prog[n_blocchi][100];
	double r_N2[n_blocchi][n_steps]; //righe=blocchi, colonne=steps, matrice in cui inserisco le medie per ogni blocco, per ogni i
	double r_N2_2[n_blocchi][n_steps]; //matrice dei quadrati dell'altra matrice

	vettore_r_N2 = new double[n_steps];
	vettore_radici = new double[n_steps];
	vettore_errori = new double[n_steps];
	x = new double[n_steps];
	y = new double[n_steps];
	z = new double[n_steps];
	n = new int[n_steps];
	
	//Parto dall'origine:
	x[0]=0;
	y[0]=0;	
	z[0]=0;

	for(int i=0; i<n_steps; i++){
		vettore_r_N2[i]=0;	
		vettore_radici[i]=0;
		vettore_errori[i]=0;
		for(int j=0; j<n_blocchi; j++){
			r_N2[j][i]=0;
			r_N2_2[j][i]=0;
			sum_prog[j][i]=0;
			sum2_prog[j][i]=0;
			err_prog[j][i]=0;
		}
	}


for (int j=0; j<n_blocchi;j++){

	for (int l=0; l<L; l++){ //ciclo sugli L RWs all'interno di ogni blocco

		for(int i=1; i<n_steps; i++){ //ciclo sugli i steps di ogni RW
			//scelgo casualmente la direzione nel continuo in cui mi devo muovere a ogni passo

			double r=0.;
     			double s=0.;
			//phi tra [0,2pi]
     			do{
       				r=rnd.Rannyu(-1,+1);
       				s=rnd.Rannyu(-1,+1);
     			}while(r*r+s*s>1 || (r==0 && s==0));//Genero due punti nella circonferenza di raggio r=1
	
     			double coseno = r/sqrt(r*r+s*s);
			double phi=0.;
			if(s<0){
				phi = 2*M_PI - acos(coseno);
			}
			else {
				phi = acos(coseno); 
			}
			//theta	tra [0,pi]
			r=0;
			s=0;
			coseno=0;
     			do{
       				r=rnd.Rannyu(-1,+1);
       				s=rnd.Rannyu(-1,+1);
     			}while(r*r+s*s>1 or (r==0 && s==0));//Genero due punti nella semi-circonferenza a destra  
     
     			coseno=(r/sqrt(r*r+s*s));
			double theta = acos(coseno); //acos ritorna un valore tra [0,pi]
			x[i] = x[i-1] + a*sin(theta)*cos(phi);
			y[i] = y[i-1] + a*sin(theta)*sin(phi);
			z[i] = z[i-1] + a*cos(theta);

		}
		
		for(int i=0; i<n_steps; i++){
			vettore_r_N2[i]=0;	
		}

		//dopo aver creato il singolo RW di 100 steps devo calcolare r_N a ogni passo
		for(int i=0; i<n_steps; i++){
			vettore_r_N2[i] = vettore_r_N2[i] + x[i]*x[i]+y[i]*y[i]+z[i]*z[i]; 
		}
		
		for(int i=0; i<n_steps; i++){ //sono all'interno del ciclo di l, quindi ogni volta somma L=100 RWs in ogni blocco
			r_N2[j][i] = r_N2[j][i] + vettore_r_N2[i];
		}
	} //fine del ciclo su l
	
	//Ora faccio la media in ogni elemento della matrice r_N2, basta dividere per L
	for(int i=0; i<n_steps; i++){ 
		r_N2[j][i] = (r_N2[j][i])/L;
		r_N2_2[j][i] = (r_N2[j][i])*(r_N2[j][i]);
	}
}	

for(int i=0; i<n_steps;i++){ //fisso il passo i-esimo(colonna) e all'interno di ogni passo faccio la somma progressiva sugli n_blocchi (righe)
	for (int j=0; j<n_blocchi; j++){
   		for (int k=0; k<j+1; k++){
			sum_prog[j][i] = sum_prog[j][i] + r_N2[k][i];
			sum2_prog[j][i] = sum2_prog[j][i]+ r_N2_2[k][i];
	
		}
		sum_prog[j][i] = sum_prog[j][i]/(j+1); //Cumulative average
    		sum2_prog[j][i] = sum2_prog[j][i]/(j+1); //Cumulative square average
   		err_prog[j][i] = error(sum_prog,sum2_prog,j,i); //Statistical uncertainty
	}
}

//La media finale per ogni i è contenuta in sum_prog[j=99][i] e il suo errore è in err_prog[j=99][i]. Ora devo fare la radice della media e trovare la sua incertezza da err_prog usando la propagazione degli errori. 

	for(int i=0; i<n_steps; i++){
		int jj=99;
		vettore_radici[i] = pow(sum_prog[jj][i],0.5);
		vettore_errori[i] = (0.5 * err_prog[jj][i] * vettore_radici[i])/(sum_prog[jj][i]);
	}	
	
	rnd.SaveSeed();

	for(int i=0; i<n_steps; i++){
		n[i]=i;
	}

	file_out.open("vettore_rN.out");
		for(int i=0; i<n_steps; i++){
			file_out << n[i] << " " << vettore_radici[i] << " " << vettore_errori[i] << endl;
		}

	file_out.close();

	return 0;
}

double error(double medie[][100], double medie2[][100], int j, int i){ //Function for statistical uncertainty estimation

	if (j==0){
	        return 0;
	}
	else{
        	return sqrt((medie2[j][i] - pow(medie[j][i],2))/(j-1));	
	 }
}
