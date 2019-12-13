#include "header.h"

static  fftw_complex   *Psi, *S, *T1, *T2, *T3;

static  double         a, b;

static  int N0;
static  int np;
static  int myid;

/* ---------------------------------------------------------------- */

void rhs_init(geom_ptr geom, phys_ptr phys,  fftw_complex *psi)
{

  int Ngrids, N;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
 
  Ngrids =  geom->Ngrids;
  N0 =  geom->N0;

  N = N0*pow(2, Ngrids-1);
  N = N*N/np;
    
  T1 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));
  T2 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));
  T3 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));

  Psi = psi;
  
  a = phys->coefA;
  b = phys->coefB;

}

/* ---------------------------------------------------------------- */


void rhs_init_S(fftw_complex *Sin)
{

  S = Sin;

}



/* ---------------------------------------------------------------- */
/* RHS is computed in place: array S holds                          */          
/* (Eta, Psi) on input and (dEta/dt, dPsi/dt) on output             */

void rhs_compute(int grid)
{

  int i, N;
 
  N = N0 * pow(2,grid);
  N = N*N/np;

  /*-- work arrays:   T1 = psi,  T2 = psi,  T3 = eta --*/

  for (i=0; i<N; i++){
    T1[i][0]  = S[i][1];
    T1[i][1]  = 0;
    T2[i][0]  = S[i][1];
    T2[i][1]  = 0;
    T3[i][0]  = S[i][0];
    T3[i][1]  = 0;
  }

  /*-- compute:   T1 = psi_x,  T2 = psi_y,  T3 = laplacian(eta) --*/
  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);
  fft_laplacian(T3, grid);

  /*--  compute dPsi = laplacian(eta) - (b/2) ( grad(psi) )^2 --*/ 
  /*--  compute [T1,T2] = (1 + a eta) grad(psi) --*/

  for (i=0; i<N; i++){

    S[i][1] = T3[i][0] - 0.5 * b * (T1[i][0]*T1[i][0] + T2[i][0]*T2[i][0]);

    T1[i][0]  = (1 + a * S[i][0]) * T1[i][0];
    T1[i][1]  = 0;
    
    T2[i][0]  = (1 + a * S[i][0]) * T2[i][0];
    T2[i][1]  = 0;
  }

  /*--  compute dEta = - div( (1 + a eta) grad(psi) )  --*/

  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);

  for (i=0; i<N; i++){

    S[i][0] = - T1[i][0] - T2[i][0];

  }

  
}


/* ---------------------------------------------------------------- */
void rhs_hamiltonian(int grid, double *Ek, double *Ep)
{
 
  int i, N, Ntot;
  
  double u,v;
  double ekin = 0;
  double epot = 0;

  
  N    = N0 * pow(2,grid);
  Ntot = N*N;
  N    = Ntot/np;

  /*-- work arrays:   T1 = eta + i psi,  T2 = eta + i psi --*/

  for (i=0; i<N; i++){
    T1[i][0]  = Psi[i][0];
    T1[i][1]  = Psi[i][1];
    T2[i][0]  = Psi[i][0];
    T2[i][1]  = Psi[i][1];
  }

  /*-- compute:   T1 = eta_x + i psi_x,  T2 = eta_y + i psi_y--*/
  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);

  
 /*-- potential energy --*/

  for (i=0; i<N; i++){
    
    u  = T1[i][0];
    v  = T2[i][0];
    epot += u*u +v*v;
    
    u  = T1[i][1];
    v  = T2[i][1];
    ekin += (u*u + v*v) * (1 + Psi[i][0]);

  }

  
  /*-- global sums --*/

  MPI_Reduce(&ekin,    &u,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&epot,    &v,       1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  *Ek = 0.5*u/Ntot;
  *Ep = 0.5*v/Ntot;

}

/* ---------------------------------------------------------------- */




