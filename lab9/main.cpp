/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <string>
#include "random.h"
#include "TSP.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes, questo" << endl;
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
   } else cerr << "PROBLEM: Unable to open seed.in, questo" << endl;

	ofstream best_path_circ, best_half_circ, final_cities_circ;
	ofstream best_path_quad, best_half_quad, final_cities_quad;

	best_path_circ.open("L1_best_path_circ.out");
	best_half_circ.open("L1_best_half_circ.out");
	final_cities_circ.open("final_cities_circ.out");
	
	best_path_quad.open("L1_best_path_quad.out");
	best_half_quad.open("L1_best_half_quad.out");
	final_cities_quad.open("final_cities_quad.out");

//CIRCONFERENZA
//Creo casualmente le 32 città sulla circonferenza.
	int N = 32;
	int Npop = 100;
	double theta = rnd.Rannyu(0,2*M_PI);
	double s = sin(theta);
	double r = cos(theta);

//Prima città, la fisso
	City inizio(r,s);
	cout << inizio.getx() << " " <<  inizio.gety() << endl;
	
//Le altre 31 città
	std::vector<City> cities;
	for (int i=0; i<N-1; i++){
		theta = rnd.Rannyu(0,2*M_PI);
		s = sin(theta);
		r = cos(theta);
		City c(r,s);
		cout << c.getx() << " " << c.gety() << endl;
		cities.push_back(c);
	}
//Creo il primo individuo
	Individuo salesman(N, inizio, cities.begin(), cities.end());
	if (salesman.check()==false){ //false per me vuol dire che non passa nella stessa città
		//cout <<"ok" << endl;
	}
	else{
		cout <<"Il percorso passa due volte nella stessa città! Non va bene!" << endl;
	}

//Creo la prima popolazione e la ordino in base a L1
	Popolazione prima(Npop, rnd, salesman); 
	prima.sort_L1();

//Generazioni
	int n_generazioni = 400;	
	for(int g=0; g<n_generazioni; g++){
		prima.new_generazione(rnd);
		//cout << "**** " << g << endl;
		prima.sort_L1();
		best_path_circ << g << " " << prima.get_pop()[0].L1() << endl;
		best_half_circ << g << " " << prima.mean_L1() << endl;
	}
	
	final_cities_circ << prima.get_pop()[0].get_partenza().getx() << " " << prima.get_pop()[0].get_partenza().gety() << endl;

	for (auto j : prima.get_pop()[0].get_percorso()){
		final_cities_circ << j.getx() << " " << j.gety() << endl;
	}
		
	cout << endl;

	for (auto i : prima.get_pop()){ 			
		cout << i.L1() << endl;
		}
	
	best_path_circ.close();
	best_half_circ.close();
	final_cities_circ.close();


	rnd.SaveSeed();
//QUADRATO

//Prima città, la fisso
	double x=rnd.Rannyu(-1,+1);
     	double y=rnd.Rannyu(-1,+1);
     	
	City inizio_q(x,y);
	cout << inizio_q.getx() << " " <<  inizio_q.gety() << endl;

//Le altre 31 città
	std::vector<City> cities_q;
	for (int i=0; i<N-1; i++){
		x=rnd.Rannyu(-1,+1);
       	y=rnd.Rannyu(-1,+1);
		City c_q(x,y);
		cout << c_q.getx() << " " << c_q.gety() << endl;
		cities_q.push_back(c_q);
	}
	
//Primo individuo
	Individuo salesman_q(N, inizio_q, cities_q.begin(), cities_q.end());
	if (salesman_q.check()==false){ //false per me vuol dire che non passa nella stessa città
		//cout <<"ok" << endl;
	}
	else{
		cout <<"Il percorso passa due volte nella stessa città! Non va bene!" << endl;
	}
	
//Prima popolazione
	Popolazione prima_q(Npop, rnd, salesman_q);
	prima_q.sort_L1();
	for (auto i : prima_q.get_pop()){ 
		cout << i.L1() << endl;
	}
	
	n_generazioni = 400;
	
	for(int g=0; g<n_generazioni; g++){
		prima_q.new_generazione(rnd);
		//cout << "**** " << g << endl;
		prima_q.sort_L1();
		best_path_quad << g << " " << prima_q.get_pop()[0].L1() << endl;
		best_half_quad << g << " " << prima_q.mean_L1() << endl;
	}
	
	final_cities_quad << prima_q.get_pop()[0].get_partenza().getx() << " " << prima_q.get_pop()[0].get_partenza().gety() << endl;

	for (auto j : prima_q.get_pop()[0].get_percorso()){
		final_cities_quad << j.getx() << " " << j.gety() << endl;
	}
		
	cout << endl;

	for (auto i : prima_q.get_pop()){ 			
		cout << i.L1() << endl;
		}
		
	best_path_quad.close();
	best_half_quad.close();
	final_cities_quad.close();
/*



	int n_generazioni = 300;
	
//Implemento la migrazione all'interno del ciclo
	int N_migr = 20;
	MPI_Status stat0_x, stat1_x, stat2_x, stat3_x, stat0_y, stat1_y, stat2_y, stat3_y;
	MPI_Request req0_x, req1_x, req2_x, req3_x, req0_y, req1_y, req2_y, req3_y;
	int itag0_x=0;
	int itag1_x=1; 
	int itag2_x=2;
	int itag3_x=3;
	int itag0_y=0;
	int itag1_y=1; 
	int itag2_y=2;
	int itag3_y=3;

	double * imesg0_x = new double[N];
	double * imesg0_y = new double[N];
	double * imesg1_x = new double[N];
	double * imesg1_y = new double[N];
	double * imesg2_x = new double[N];
	double * imesg2_y = new double[N];
	double * imesg3_x = new double[N];
	double * imesg3_y = new double[N];

	int i=1;
	for(int g=0; g<n_generazioni; g++){
	
	MPI_Barrier(MPI_COMM_WORLD); 
		if((g+1)%N_migr==0){
			cout << "Migrazione numero: " << (g+1)/N_migr << endl;
			int index_send;
			int index_altro;
			int index_altro1;
			if(rank==0){
				imesg0_x[0] = prima_q.get_pop()[0].get_partenza().getx();
				imesg0_y[0] = prima_q.get_pop()[0].get_partenza().gety();
				i=1;
				for (auto j : prima_q.get_pop()[0].get_percorso()){
					imesg0_x[i] =  j.getx();	
					imesg0_y[i] =  j.gety();	
					i=i+1;	
				}
			} else if(rank==1){
				imesg1_x[0] = prima_q.get_pop()[0].get_partenza().getx();
				imesg1_y[0] = prima_q.get_pop()[0].get_partenza().gety();
				i=1;
				for (auto j : prima_q.get_pop()[0].get_percorso()){
					imesg1_x[i] =  j.getx();	
					imesg1_y[i] =  j.gety();	
					i=i+1;
				}
			} else if(rank==2){
				imesg2_x[0] = prima_q.get_pop()[0].get_partenza().getx();
				imesg2_y[0] = prima_q.get_pop()[0].get_partenza().gety();
				i=1;
				for (auto j : prima_q.get_pop()[0].get_percorso()){
					imesg2_x[i] =  j.getx();	
					imesg2_y[i] =  j.gety();	
					i=i+1;
				}
			} else if(rank==3){
				imesg3_x[0] = prima_q.get_pop()[0].get_partenza().getx();
				imesg3_y[0] = prima_q.get_pop()[0].get_partenza().gety();
				i=1;
				for (auto j : prima_q.get_pop()[0].get_percorso()){
					imesg3_x[i] =  j.getx();	
					imesg3_y[i] =  j.gety();	
					i=i+1;
				}
			}

			do{
				index_send=int(rnd.Rannyu(0,4));
			} while(index_send==0);

			do{
				index_altro=int(rnd.Rannyu(0,4));
			} while(index_altro==0 || index_altro==index_send);

			do{
				index_altro1=int(rnd.Rannyu(0,4));
			} while(index_altro1==0 || index_altro1==index_send || index_altro1==index_altro);
			
			cout << " index" << index_send << endl;
			cout << " index altro " << index_altro << endl;
			cout << " index altro1 " << index_altro1 << endl;

			if(rank==0){
				MPI_Isend(&imesg0_x[0],N,MPI_DOUBLE_PRECISION,index_send,itag0_x,MPI_COMM_WORLD,&req0_x);
				MPI_Isend(&imesg0_y[0],N,MPI_DOUBLE_PRECISION,index_send,itag0_y,MPI_COMM_WORLD,&req0_y);
				cout << " index" << index_send << endl;
				if(index_send==1){
					MPI_Recv(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_send,itag1_x, MPI_COMM_WORLD,&stat1_x);
					MPI_Recv(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_send,itag1_y, MPI_COMM_WORLD,&stat1_y);
					prima_q.migrazione(imesg1_x,imesg1_y);
					prima_q.sort_L1();
				}else if(index_send==2){
					MPI_Recv(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_send,itag2_x, MPI_COMM_WORLD,&stat2_x);
					MPI_Recv(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_send,itag2_y, MPI_COMM_WORLD,&stat2_y);
					prima_q.migrazione(imesg2_x,imesg2_y);
					prima_q.sort_L1();
				}else if(index_send==3){
					MPI_Recv(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_send,itag3_x, MPI_COMM_WORLD,&stat3_x);
					MPI_Recv(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_send,itag3_y, MPI_COMM_WORLD,&stat3_y);
					prima_q.migrazione(imesg3_x,imesg3_y);
					prima_q.sort_L1();
				}
				//cout<<"messaggio ricevuto da 0 = "<<imesg2_x[10]<<endl;}
				
			}else if(rank==index_send){
				if(rank==1){
					MPI_Isend(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,0,itag1_x, MPI_COMM_WORLD,&req1_x);
					MPI_Isend(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,0,itag1_y, MPI_COMM_WORLD,&req1_y);
					MPI_Recv(&imesg0_x[0],N,MPI_DOUBLE_PRECISION,0,itag0_x, MPI_COMM_WORLD,&stat0_x);
					MPI_Recv(&imesg0_y[0],N,MPI_DOUBLE_PRECISION,0,itag0_y, MPI_COMM_WORLD,&stat0_y);
				}else if(rank==2){
					MPI_Isend(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,0,itag2_x, MPI_COMM_WORLD,&req2_x);
					MPI_Isend(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,0,itag2_y, MPI_COMM_WORLD,&req2_y);
					MPI_Recv(&imesg0_x[0],N,MPI_DOUBLE_PRECISION,0,itag0_x, MPI_COMM_WORLD,&stat0_x);
					MPI_Recv(&imesg0_y[0],N,MPI_DOUBLE_PRECISION,0,itag0_y, MPI_COMM_WORLD,&stat0_y);
				}else if(rank==3){
					MPI_Isend(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,0,itag3_x, MPI_COMM_WORLD,&req3_x);
					MPI_Isend(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,0,itag3_y, MPI_COMM_WORLD,&req3_y);
					MPI_Recv(&imesg0_x[0],N,MPI_DOUBLE_PRECISION,0,itag0_x, MPI_COMM_WORLD,&stat0_x);
					MPI_Recv(&imesg0_y[0],N,MPI_DOUBLE_PRECISION,0,itag0_y, MPI_COMM_WORLD,&stat0_y);
				}
				prima_q.migrazione(imesg0_x,imesg0_y);
				prima_q.sort_L1();
				//cout<<"messaggio = "<<imesg[0]<<endl;}
			}	

			if(rank==index_altro){
				if(rank==1){
					MPI_Isend(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag1_x, MPI_COMM_WORLD,&req1_x);
					MPI_Isend(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag1_y, MPI_COMM_WORLD,&req1_y);
					if(index_altro1==2){
						MPI_Recv(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag2_x, MPI_COMM_WORLD,&stat2_x);
						MPI_Recv(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag2_y, MPI_COMM_WORLD,&stat2_y);
						prima_q.migrazione(imesg2_x,imesg2_y);
						prima_q.sort_L1();
					}
					else if(index_altro1==3){
						MPI_Recv(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag3_x, MPI_COMM_WORLD,&stat3_x);
						MPI_Recv(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag3_y, MPI_COMM_WORLD,&stat3_y);
						prima_q.migrazione(imesg3_x,imesg3_y);
						prima_q.sort_L1();
					}
				}else if(rank==2){
					MPI_Isend(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag2_x, MPI_COMM_WORLD,&req2_x);
					MPI_Isend(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag2_y, MPI_COMM_WORLD,&req2_y);
					if(index_altro1==1){
						MPI_Recv(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag1_x, MPI_COMM_WORLD,&stat1_x);
						MPI_Recv(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag1_y, MPI_COMM_WORLD,&stat1_y);
						prima_q.migrazione(imesg1_x,imesg1_y);
						prima_q.sort_L1();
					}
					else if(index_altro1==3){
						MPI_Recv(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag3_x, MPI_COMM_WORLD,&stat3_x);
						MPI_Recv(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag3_y, MPI_COMM_WORLD,&stat3_y);
						prima_q.migrazione(imesg3_x,imesg3_y);
						prima_q.sort_L1();
					}
				}else if(rank==3){
					MPI_Isend(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag3_x, MPI_COMM_WORLD,&req3_x);
					MPI_Isend(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag3_y, MPI_COMM_WORLD,&req3_y);
					if(index_altro1==1){
						MPI_Recv(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag1_x, MPI_COMM_WORLD,&stat1_x);
						MPI_Recv(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag1_y, MPI_COMM_WORLD,&stat1_y);
						prima_q.migrazione(imesg1_x,imesg1_y);
						prima_q.sort_L1();
					}
					else if(index_altro1==2){
						MPI_Recv(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag2_x, MPI_COMM_WORLD,&stat2_x);
						MPI_Recv(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_altro1,itag2_y, MPI_COMM_WORLD,&stat2_y);
						prima_q.migrazione(imesg2_x,imesg2_y);
						prima_q.sort_L1();
					}
				}
				//cout<<"messaggio = "<<imesg[0]<<endl;}	
			}else if(rank==index_altro1){
				if(rank==1){
					MPI_Isend(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag1_x, MPI_COMM_WORLD,&req1_x);
					MPI_Isend(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag1_y, MPI_COMM_WORLD,&req1_y);
					if(index_altro==2){
						MPI_Recv(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag2_x, MPI_COMM_WORLD,&stat2_x);
						MPI_Recv(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag2_y, MPI_COMM_WORLD,&stat2_y);
						prima_q.migrazione(imesg2_x,imesg2_y);
						prima_q.sort_L1();
					}
					else if(index_altro==3){
						MPI_Recv(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag3_x, MPI_COMM_WORLD,&stat3_x);
						MPI_Recv(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag3_y, MPI_COMM_WORLD,&stat3_y);
						prima_q.migrazione(imesg3_x,imesg3_y);
						prima_q.sort_L1();
					}
				}else if(rank==2){
					MPI_Isend(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag2_x, MPI_COMM_WORLD,&req2_x);
					MPI_Isend(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag2_y, MPI_COMM_WORLD,&req2_y);
					if(index_altro==1){
						MPI_Recv(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag1_x, MPI_COMM_WORLD,&stat1_x);
						MPI_Recv(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag1_y, MPI_COMM_WORLD,&stat1_y);
						prima_q.migrazione(imesg1_x,imesg1_y);
						prima_q.sort_L1();
					}
					else if(index_altro==3){
						MPI_Recv(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag3_x, MPI_COMM_WORLD,&stat3_x);
						MPI_Recv(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag3_y, MPI_COMM_WORLD,&stat3_y);
						prima_q.migrazione(imesg3_x,imesg3_y);
						prima_q.sort_L1();
					}
				}else if(rank==3){
					MPI_Isend(&imesg3_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag3_x, MPI_COMM_WORLD,&req3_x);
					MPI_Isend(&imesg3_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag3_y, MPI_COMM_WORLD,&req3_y);
					if(index_altro==1){
						MPI_Recv(&imesg1_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag1_x, MPI_COMM_WORLD,&stat1_x);
						MPI_Recv(&imesg1_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag1_y, MPI_COMM_WORLD,&stat1_y);
						prima_q.migrazione(imesg1_x,imesg1_y);
						prima_q.sort_L1();
					}
					else if(index_altro==2){
						MPI_Recv(&imesg2_x[0],N,MPI_DOUBLE_PRECISION,index_altro,itag2_x, MPI_COMM_WORLD,&stat2_x);
						MPI_Recv(&imesg2_y[0],N,MPI_DOUBLE_PRECISION,index_altro,itag2_y, MPI_COMM_WORLD,&stat2_y);
						prima_q.migrazione(imesg2_x,imesg2_y);
						prima_q.sort_L1();
					}
				}
			}

		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		prima_q.new_generazione(rnd_new);
		//cout << "**** " << g << endl;
		prima_q.sort_L1();
		if(rank==0){
			best_path_quad0 << g << " " << prima_q.get_pop()[0].L1() << endl;
			best_half_quad0 << g << " " << prima_q.mean_L1() << endl;
		} else if(rank==1){
			best_path_quad1 << g << " " << prima_q.get_pop()[0].L1() << endl;
			best_half_quad1 << g << " " << prima_q.mean_L1() << endl;
		} else if(rank==2){
			best_path_quad2 << g << " " << prima_q.get_pop()[0].L1() << endl;
			best_half_quad2 << g << " " << prima_q.mean_L1() << endl;
		} else if(rank==3){
			best_path_quad3 << g << " " << prima_q.get_pop()[0].L1() << endl;
			best_half_quad3 << g << " " << prima_q.mean_L1() << endl;
		}
		//for (auto i : prima.get_pop()){ 
		//	cout << i.L1() << endl;
		//}
		
	}
	
	if(rank==0){
		final_cities_quad0 << prima_q.get_pop()[0].get_partenza().getx() << " " << prima_q.get_pop()[0].get_partenza().gety() << endl;
	} else if(rank==1){
		final_cities_quad1 << prima_q.get_pop()[0].get_partenza().getx() << " " << prima_q.get_pop()[0].get_partenza().gety() << endl;
	} else if(rank==2){
		final_cities_quad2 << prima_q.get_pop()[0].get_partenza().getx() << " " << prima_q.get_pop()[0].get_partenza().gety() << endl;
	} else if(rank==3){
		final_cities_quad3 << prima_q.get_pop()[0].get_partenza().getx() << " " << prima_q.get_pop()[0].get_partenza().gety() << endl;
	}

	if(rank==0){
		for (auto j : prima_q.get_pop()[0].get_percorso()){
			final_cities_quad0 << j.getx() << " " << j.gety() << endl;
		}
	} else if(rank==1){
		for (auto j : prima_q.get_pop()[0].get_percorso()){
			final_cities_quad1 << j.getx() << " " << j.gety() << endl;
		}
	} else if(rank==2){
		for (auto j : prima_q.get_pop()[0].get_percorso()){
			final_cities_quad2 << j.getx() << " " << j.gety() << endl;
		}
	} else if(rank==3){		
		for (auto j : prima_q.get_pop()[0].get_percorso()){
			final_cities_quad3 << j.getx() << " " << j.gety() << endl;
		}
	}		
	cout << endl;
	if(rank==0){
		for (auto i : prima_q.get_pop()){ 
			cout << i.L1() << endl;
		}
	}
	best_path_quad0.close();
	best_half_quad0.close();
	final_cities_quad0.close();
	best_path_quad1.close();
	best_half_quad1.close();
	final_cities_quad1.close();
	best_path_quad2.close();
	best_half_quad2.close();
	final_cities_quad2.close();
	best_path_quad3.close();
	best_half_quad3.close();
	final_cities_quad3.close();
	rnd.SaveSeed();
	rnd_new.SaveSeed();

	MPI_Finalize();
*/
   	return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
