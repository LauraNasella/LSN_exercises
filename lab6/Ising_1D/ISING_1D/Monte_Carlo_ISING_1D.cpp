/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  ifstream ReadInput;
  ReadInput.open("input.dat");
  ReadInput >> restart;
  ReadInput.close();
  if (restart == 0){ //prima volta in cui faccio girare il codice
  	double tem=2.0;
  	Input(tem,restart); //Inizialization
  	for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  	{
    		Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep)
    		{
     			Move(metro);
     			Measure(iblk,restart);
      			Accumulate(); //Update block averages
    		}
    		Averages(iblk);   //Print results for current block
  	}
  ConfFinal(); //Write final configuration
  }
  else{
  	for(int t=0; t<31; t++){
  		double tem = 2.0;
  		tem = tem - t*0.05;
  		Input(tem,restart); //Inizialization
  		for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  		{
    			Reset(iblk);   //Reset block averages
    			for(int istep=1; istep <= nstep; ++istep)
    			{
     				Move(metro);
     				Measure(iblk,restart);
      				Accumulate(); //Update block averages
    			}
    			Averages(iblk);   //Print results for current block
  		}
  	ConfFinal(); //Write final configuration
  	}
  }
  return 0;
}

void Input(double tem,int restart) //Ho modificato input in modo che prenda un double, cioè la temperatura
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  temp=tem;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> restart; //Variabile aggiunta che assume come valori 0 (primo giro del codice) o 1 (permette la ripartenza successiva)

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  
  ReadInput.close();

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

	if(restart == 0){ //Se è 0 significa che sto facendo girare il codice per la prima volta, quindi prepara la configurazione iniziale degli spin in modo casuale
		//initial configuration
  		for (int i=0; i<nspin; ++i){
    			if(rnd.Rannyu() >= 0.5) s[i] = 1;
    			else s[i] = -1;
  		}
  	}
  	else{ //Cioè riparto dalla configurazione precedente di spin, salvata in config.final
  		ifstream ReadConfig;
  		ReadConfig.open("config.final");
  		for (int i=0; i<nspin; ++i){
    			ReadConfig >> s[i];
  		}
  		ReadConfig.close();
  	}
  	  		
  	//Evaluate energy etc. of the initial configuration
  	int iblk=0;
  	Measure(iblk,restart);
  		
  	//Print initial values for the potential energy and virial
  	cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
  int o;
  //double p;
  double energy_old, energy_new, sm, deltaE;
  //double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    	if(metro==1) //Metropolis
    	{
// INCLUDE YOUR CODE HERE
		sm = s[o]*(-1); //cioè preso lo spin o in modo casuale, propongo la mossa di flipparlo
		energy_old = Boltzmann(s[o],o);
		energy_new = Boltzmann(sm,o);
		deltaE = energy_new - energy_old;
	
		double alpha = 0.;
	
		if(deltaE < 0){
			alpha = 1;
		}
		else{
			alpha = exp(-deltaE*beta);
		}
	
		double r = rnd.Rannyu(); //genero un numero casuale uniformemente tra 0 e 1
		if(r<=alpha){
			attempted = attempted +1;
   			accepted = accepted +1;
			s[o]=sm;
		}	
		else{
			attempted = attempted +1;
		}		
	}

  	else //Gibbs sampling
    	{
// INCLUDE YOUR CODE HERE
		deltaE = -2.*Boltzmann(1,o); 
		double alpha = 1./(1.+exp(-beta*deltaE));
		double r = rnd.Rannyu(); //genero un numero casuale uniformemente tra 0 e 1
		
	 	if(r<=alpha){
	    		s[o]=+1;
	  	}
	  	else{
	    		s[o]=-1;
		}
	}
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure(int iblk,int restart)
{
  //int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
// INCLUDE YOUR CODE HERE
	m += s[i];
	//cout << u << endl;
  }
  walker[iu] = u;
  //cout << walker[iu] << endl;
// INCLUDE YOUR CODE HERE
	walker[ic] = u*u; //mi serve per la capacità termica, che poi calcolo effettivamente in Averages
	walker[im] = m;
	walker[ix] = m*m;
	
	const int wd=12;
	ofstream Mag_eq;
	if(h==0.02){
    		if(temp==0.5){ 
    			if(iblk==1){
    				if(metro==1){
    					Mag_eq.open("equil.mag.metro",ios::app);
    				}
    				else{
    					Mag_eq.open("equil.mag.gibbs",ios::app);
    				}
    				Mag_eq << setw(wd) << iblk << setw(wd)<< temp << setw(wd) << walker[im] << endl;  
    			}
    		}
    	}
    	Mag_eq.close();
    	
    	ofstream Ene_eq;
    	if(restart == 0){
	if(h==0.0){
    		if(temp==2){ 
    			if(iblk==1){
    				if(metro==1){
    					Ene_eq.open("equil.ene.metro",ios::app);
    				}
    				else{
    					Ene_eq.open("equil.ene.gibbs",ios::app);
    				}
    				Ene_eq << setw(wd) << iblk << setw(wd)<< temp << setw(wd) << walker[iu] << endl;  
    			}
    		}
    	}
    	}
    	Ene_eq.close();
    	
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    if(metro==1){
    	Ene.open("output.ene.metro.0",ios::app);
    }
    else{
    	Ene.open("output.ene.gibbs.0",ios::app);
    }
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy, già divisa per N
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    if(h==0){
   	if(iblk == nblk){
		Ene << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	}
    }
    Ene.close();

// INCLUDE YOUR CODE HERE

	if(metro==1){
    		Heat.open("output.heat.metro.0",ios::app);
    	}
    	else{
    		Heat.open("output.heat.gibbs.0",ios::app);
    	}
    
	stima_c = (pow(beta,2)*((blk_av[ic]/blk_norm)-pow(blk_av[iu]/blk_norm,2)))/(double)nspin; //(blk_av[ic]/blk_norm) è <H^2>, trovo C già divisa per N, cioè il calore specifico
	glob_av[ic] += stima_c;
	glob_av2[ic] += stima_c*stima_c;
    	err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    	if(h==0){
    		if(iblk == nblk){
    			Heat << setw(wd) << iblk <<  setw(wd) << temp << setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    		}
    	}
    	Heat.close();
    	
    	if(metro==1){
    		Mag.open("output.mag.metro.0",ios::app);
    	}
    	else{
    		Mag.open("output.mag.gibbs.0",ios::app);
    	}
    	
	stima_m = blk_av[im]/blk_norm/(double)nspin; //già divisa per N
	glob_av[im] += stima_m;
	glob_av2[im] += stima_m*stima_m;
    	err_m=Error(glob_av[im],glob_av2[im],iblk);
    	if(h==0.02){
    		if(iblk == nblk){
    			Mag << setw(wd) << iblk << setw(wd) << temp  << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    		}
    	}
    	Mag.close();
    	
    	//Queste formule valgono solo per h=0, cioè pongo direttamente (il valor medio della magnetizzazione)^2 = 0.
    	if (h==0){
    		if(metro==1){
    			Chi.open("output.chi.metro.0",ios::app);
    		}
    		else{
    			Chi.open("output.chi.gibbs.0",ios::app);
    		}
		stima_x = beta*blk_av[ix]/blk_norm/(double)nspin;
		glob_av[ix] += stima_x;
		glob_av2[ix] += stima_x*stima_x;
    		err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    		if(iblk == nblk){
    			Chi << setw(wd) << iblk << setw(wd)<< temp << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    		}
    		Chi.close();
    	}
    	
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
