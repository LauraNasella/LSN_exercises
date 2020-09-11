#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>

using namespace std;

double stat_unc(double * istant, int N, int M);

int main(){

	ifstream pres_ist, epot_ist;
	pres_ist.open("output_pres.dat");
  	epot_ist.open("output_epot.dat");
  	
  	ofstream fileout;
  	fileout.open("inc_stat.dat",ios::app);
  	
  	int M=500000;
  	int L_min = 10;
  	int L_max = 5000; 
  	int N_min = int(M/L_max); //numero di blocchi
  	int N_max = int(M/L_min);
  	
	double * U = new double[M];
	double * P = new double[M];
	
	for (int i=0; i<M; i++){
		pres_ist >> P[i];
		epot_ist >> U[i];
	}
	
	pres_ist.close();
	epot_ist.close();
        
        int conta = 0;
        
	for (int N=N_min; N<N_max+1; N++){
		if((M%N)==0){	
			conta++;
			cout << N << " " << conta << endl;
			double err_U = stat_unc(U,N,M);
			double err_P = stat_unc(P,N,M);
		
        		fileout << int(M/N) << " " << err_U << " " << err_P << endl;
        	}
	}

	fileout.close();
	delete[] U;
	delete[] P;
	
	return 0;
}

double stat_unc(double * istant, int N, int M){
    
	int L=int(M/N);

	double * mean = NULL;
	double * mean2 = NULL;
	mean = new double[N];
	mean2 = new double[N];
       	 	
        for(int i=0; i<N; i++){
        	double sum = 0;
        	for(int j=0; j<L; j++){
        		int k = j+i*L;
			sum = sum + istant[k];
		}
		mean[i] = sum/L;
		mean2[i] = mean[i]*mean[i];
	}
	
	double media = 0;
	double media2 = 0;
	
	for(int i=0; i<N; i++){
		media = media + mean[i];
		media2 = media2 + mean2[i];
	}
	media = media/N;
	media2 = media2/N;

	return sqrt((media2 - media*media)/(N-1));

}
