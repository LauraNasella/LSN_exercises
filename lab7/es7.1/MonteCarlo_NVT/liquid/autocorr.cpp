#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

using namespace std;

double * autocorr(double * ist,int steps, int M);

int main(){

	ifstream pres_ist, epot_ist;
	pres_ist.open("output_pres.dat");
  	epot_ist.open("output_epot.dat");
  	
  	int M=500000;
  	
	double * U = new double[M];
	double * P = new double[M];
	
	for (int i=0; i<M; i++){
		pres_ist >> P[i];
		epot_ist >> U[i];
	}
	
	pres_ist.close();
	epot_ist.close();

	int steps = 400;
	double * AC_pres = autocorr(P,steps,M);
	double * AC_epot = autocorr(U,steps,M);
	
	ofstream Autoc_pres, Autoc_epot;
	Autoc_pres.open("autoc_pres.dat");
	Autoc_epot.open("autoc_epot.dat");
	for(int i=0; i<steps; i++){
		Autoc_pres << AC_pres[i] << endl;
		Autoc_epot << AC_epot[i] << endl;
	}
	
	return 0;
}

double * autocorr(double * ist, int steps, int M){

	double * AC = new double[steps];

	double varianza_ist = 0;
	double * ist2 = new double[M];
	double media = 0;
	double media2 = 0;
	
	for(int i=0; i<M;i++){
		ist2[i]= ist[i]*ist[i];
	}
	for(int i=0; i<M; i++){
		media = media + ist[i];
		media2 = media2 + ist2[i];
	}

	varianza_ist = media2/M - pow((media/M),2);
	
	for(int tau=0; tau<steps; tau++){
		double mean_t = 0;
		double mean_ttau = 0;
		double mean_prodotto = 0;
		for(int t=0; t<(M-tau); t++){
			mean_t = mean_t + ist[t];
			mean_ttau = mean_ttau + ist[t+tau];
			mean_prodotto = mean_prodotto + ist[t]*ist[t+tau];
		}
		AC[tau] = (mean_prodotto/(M-tau) - (mean_t)*(mean_ttau)/pow((M-tau),2))/varianza_ist;
	}
	return AC;
}	
		
