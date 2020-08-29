#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"
#include <cmath>

using namespace std;

int main (int argc, char *argv[]){

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

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
   
	int n=10000;
	int N[4]={1,2,10,100}; 
     
	//Vettori che contengono gli S_N per i 4 casi
  	double * dado_std[4] ;
   	double * dado_exp[4] ;
   	double * dado_lor[4] ;

	for(int i=0; i<4; i++){
   		dado_std[i]=new double[n];
   		dado_exp[i]=new double[n];
   		dado_lor[i]=new double[n];
   	}
   
   	//Riempio i vettori con le 3 distribuzioni nei 4 casi 
	for(int i=0; i<4; i++){
     		for(int j=0; j<n; j++){
			dado_std[i][j]=0.;
      			dado_exp[i][j]=0.;
       			dado_lor[i][j]=0.;
      			for(int k=0; k<N[i]; k++){
	 			dado_std[i][j]= dado_std[i][j] + int(rnd.Rannyu(0,10));
	 			dado_exp[i][j]= dado_exp[i][j] + rnd.Exp(1);
				dado_lor[i][j]= dado_lor[i][j] + rnd.Lorentzian(0,1);
       			}
       				dado_std[i][j]=dado_std[i][j]/N[i];
       				dado_exp[i][j]=dado_exp[i][j]/N[i];
      				dado_lor[i][j]=dado_lor[i][j]/N[i];
     		}
   	}

   	//Dado standard
   	ofstream file_out;
   	file_out.open("dado_std.dat");
   
	for(int j=0; j<n; j++){
     		for(int i=0; i<4; i++){
			file_out << dado_std[i][j]<<" ";
		}
     		file_out<<endl;
   	}
   	file_out.close();

	//Dado esponenziale
   	file_out.open("dado_exp.dat");
   
	for(int j=0; j<n; j++){
		for(int i=0; i<4; i++){
			file_out << dado_exp[i][j]<<" ";
     		}
     		file_out<<endl;
   	}
   	file_out.close();

	//Dado Lorentziano
   	file_out.open("dado_lor.dat");
   
	for(int j=0; j<n; j++){
		for(int i=0; i<4; i++){
			file_out << dado_lor[i][j]<<" ";
     		}
     		file_out<<endl;
   	}
   	file_out.close();
      
	rnd.SaveSeed();

	return 0;
}



