/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=105;
int n_props;
int iv,ik,it,ie,ip, igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres, bin_size, nbins,err_gdir;
double ptail;
double walker[m_props];
double glob_av[m_props],glob_av2[m_props];
double blk_av[m_props],blk_norm,accepted,attempted;
const int nblk = 50;
// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double vx_2[m_part],vy_2[m_part],vz_2[m_part];

// thermodynamical state
int npart;
int restart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
void Average(int index);
double error(double * medie, double * medie2, int i);
double Error(double,double,int);
void Reset(int);
void Accumulate(void);
void Averages(int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
