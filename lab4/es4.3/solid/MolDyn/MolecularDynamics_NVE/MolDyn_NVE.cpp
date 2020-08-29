/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  ConfFinal();         //Write final configuration to restart
  Average(0);
  Average(1);
  Average(2);
  Average(3);
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> restart; //Ho aggiunto questa variabile che assume come valori 0 (primo giro del codice) o 1 (permette la ripartenza successiva)
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

	if(restart == 0){ //Se è 0 significa che sto facendo girare il codice per la prima volta, quindi prepara le velocità iniziali in modo casuale
	//Read initial configuration
  		cout << "Read initial configuration from file config.0 " << endl << endl;
  		ReadConf.open("config.0");
  		for (int i=0; i<npart; ++i){
 			ReadConf >> x[i] >> y[i] >> z[i];
    			x[i] = x[i] * box;
   			y[i] = y[i] * box;
    			z[i] = z[i] * box;
  		}
  		ReadConf.close();

	//Prepare initial velocities
 	 	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   		double sumv[3] = {0.0, 0.0, 0.0};
   		for (int i=0; i<npart; ++i){
     			vx[i] = rand()/double(RAND_MAX) - 0.5;
     			vy[i] = rand()/double(RAND_MAX) - 0.5;
     			vz[i] = rand()/double(RAND_MAX) - 0.5;

     			sumv[0] += vx[i];
     			sumv[1] += vy[i];
     			sumv[2] += vz[i];
   		}
   		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
  		double sumv2 = 0.0, fs;
   		for (int i=0; i<npart; ++i){
     			vx[i] = vx[i] - sumv[0];
     			vy[i] = vy[i] - sumv[1];
     			vz[i] = vz[i] - sumv[2];

     			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   		}
   		sumv2 /= (double)npart;

   		fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
  		for (int i=0; i<npart; ++i){
     			vx[i] *= fs;
     			vy[i] *= fs;
     			vz[i] *= fs;

     			xold[i] = Pbc(x[i] - vx[i] * delta);
     			yold[i] = Pbc(y[i] - vy[i] * delta);
     			zold[i] = Pbc(z[i] - vz[i] * delta);
   		}
   	} 
   	else{
   	//Read initial and old configuration
  		cout << "Read initial configuration from file config.0 " << endl << endl;
  		ReadConf.open("config.0");
  		for (int i=0; i<npart; ++i){
 			ReadConf >> x[i] >> y[i] >> z[i];
    			x[i] = x[i] * box;
   			y[i] = y[i] * box;
    			z[i] = z[i] * box;
  		}
  		ReadConf.close();
		cout << "Read old configuration from file old.0 " << endl << endl;
  		ReadConf.open("old.0");
  		for (int i=0; i<npart; ++i){
 			ReadConf >> xold[i] >> yold[i] >> zold[i];
    			xold[i] = xold[i] * box;
   			yold[i] = yold[i] * box;
    			zold[i] = zold[i] * box;
  		}
  		ReadConf.close();

	//Faccio uno step dell'algoritmo di Verlet per trovare xnew,ynew,znew e da questo e x,y,z trovare le velocità e la temperatura
		double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
		double sumv2 = 0.0;
		double actual_temp = 0.0;
		double fs;
		
  		for(int i=0; i<npart; ++i){ //Force acting on particle i
    			fx[i] = Force(i,0);
   			fy[i] = Force(i,1);
    			fz[i] = Force(i,2);
  		}
  		
		for(int i=0; i<npart; ++i){ //Verlet integration scheme

    			xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    			ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    			znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    			vx_2[i] = Pbc(xnew - x[i])/(delta);	
    			vy_2[i] = Pbc(ynew - y[i])/(delta);
   			vz_2[i] = Pbc(znew - z[i])/(delta);
   			
   			sumv2 += vx_2[i]*vx_2[i] + vy_2[i]*vy_2[i] + vz_2[i]*vz_2[i]; 
		}
  		sumv2 /= (double)npart;
   		actual_temp = sumv2/3;
   		
   		fs = sqrt(temp / actual_temp);   // fs = velocity scale factor 
  		for (int i=0; i<npart; ++i){
     			vx[i] = vx_2[i] * fs;
     			vy[i] = vy_2[i] * fs;
     			vz[i] = vz_2[i] * fs;

     			xold[i] = Pbc(x[i] - vx[i] * delta);
     			yold[i] = Pbc(y[i] - vy[i] * delta);
     			zold[i] = Pbc(z[i] - vz[i] * delta);
   		}
   	}
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void ConfFinal(void){ //Write final and old configuration (modificata da me)
  ofstream WriteConf;
  ofstream WriteConf1;

  cout << "Print final configuration to file config.0 and and old spatial configuration to file old.0" << endl << endl;
  WriteConf.open("config.0");
  WriteConf1.open("old.0");

  for (int i=0; i<npart; ++i){
    	WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  for (int i=0; i<npart; ++i){
    	WriteConf1 << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  WriteConf1.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

void Average(int index){
	int nblocks = 50;
	int nmeasure = nstep/10;
	int L = nmeasure/nblocks;
	
	double * prop_ist = NULL;
	double * prop = NULL;
	double * prop2 = NULL;
	double * sum_prog = NULL;
	double * sum2_prog = NULL;
	double * err_prog = NULL;
	int * x = NULL;
	
	prop_ist = new double[nmeasure];
	prop = new double[nblocks];
	prop2 = new double[nblocks];
	sum_prog = new double[nblocks];
	sum2_prog = new double[nblocks];
	err_prog = new double[nblocks];
	x = new int[nblocks];
	
	if (index == 0){ //En.Potenziale
		ifstream Epot_ist;
  		Epot_ist.open("output_epot.dat");
  		for(int i=0; i<nmeasure ; i++){
  			Epot_ist >> prop_ist[i];
  		}
  		Epot_ist.close();
  	}
  	
  	if (index == 1){ //En.Cinetica
  		ifstream Ekin_ist;
  		Ekin_ist.open("output_ekin.dat");
  		for(int i=0; i<nmeasure ; i++){
  			Ekin_ist >> prop_ist[i];
  		}
  		Ekin_ist.close();
  	}
  	
  	if (index == 2){ //En.Totale
  		ifstream Etot_ist;
  		Etot_ist.open("output_etot.dat");
  		for(int i=0; i<nmeasure ; i++){
  			Etot_ist >> prop_ist[i];
  		}
  		Etot_ist.close();
  	}
  	
    	if (index == 3){ //Temperatura
    		ifstream Temp_ist;
  		Temp_ist.open("output_temp.dat");
  		for(int i=0; i<nmeasure ; i++){
  			Temp_ist >> prop_ist[i];
  		}
  		Temp_ist.close();
  	}
  	
  	for (int i=0; i<nblocks; i++){
    		double sum = 0;
   		for (int j=0; j<L; j++){
        		int k = j+i*L;
			sum = sum + prop_ist[k];
		}
		prop[i] = sum/L;     
    		prop2[i] = pow(prop[i],2);
	}

	for (int i=0; i<nblocks; i++){
   		for (int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + prop[j];
			sum2_prog[i] = sum2_prog[i]+ prop2[j];
	
		}
		sum_prog[i] = sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i] = sum2_prog[i]/(i+1); //Cumulative square average
   		err_prog[i] = error(sum_prog,sum2_prog,i); //Statistical uncertainty
	}
	
	for(int i=0; i<nblocks; i++){
		x[i]=i;
	}
	
	if (index == 0){
		ofstream Epot_ave;
		Epot_ave.open("ave_epot.dat",ios::app);
		for(int i=0; i<nblocks; i++){
			Epot_ave << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
		}
		Epot_ave.close();
	}
	
	if (index == 1){
		ofstream Ekin_ave;
		Ekin_ave.open("ave_ekin.dat",ios::app);
		for(int i=0; i<nblocks; i++){
			Ekin_ave << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
		}
		Ekin_ave.close();
	}
	
	if (index == 2){
		ofstream Etot_ave;
	  	Etot_ave.open("ave_etot.dat",ios::app);
		for(int i=0; i<nblocks; i++){
			Etot_ave << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
		}
		Etot_ave.close();
	}
	
	if (index == 3){
		ofstream Temp_ave;
	  	Temp_ave.open("ave_temp.dat",ios::app);
		for(int i=0; i<nblocks; i++){
			Temp_ave << x[i] << " " << sum_prog[i] << " " << err_prog[i] << endl;
		}
		Temp_ave.close();
	}
	
}

double error(double * medie, double * medie2, int i){ 

	if (i==0){
	        return 0;
	}
	else{
        	return sqrt((medie2[i] - pow(medie[i],2))/(i-1));	
	 }
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
