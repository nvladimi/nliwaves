#include "header.h"

static  fftw_complex   *Psi, *T1, *T2, *T3;

static  double         a, b;

static  int N0;
static  int np;

/* ---------------------------------------------------------------- */

void rhs_init(geom_ptr geom, phys_ptr phys)
{

  int Ngrids, N;

  MPI_Comm_size(MPI_COMM_WORLD, &np);
 
  Ngrids =  geom->Ngrids;
  N0 =  geom->N0;

  N = N0*pow(2, Ngrids-1);
  N = N*N/np;
    
  T1 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));
  T2 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));
  T3 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));

  a = phys->coefA;
  b = phys->coefB;

}

/* ---------------------------------------------------------------- */


void rhs_init_S(fftw_complex *psi)
{

  Psi = psi;

}



/* ---------------------------------------------------------------- */
void rhs_compute(int grid)   // Eta, Psi = rhs(Eta, Psi) 
{

  int i, N;
 
  N = N0 * pow(2,grid);
  N = N*N/np;

  /*-- work arrays:   T1 = v, T2 = v, T3 = u --*/

  for (i=0; i<N; i++){
    T1[i][0]  = - Psi[i][1];
    T1[i][1]  = 0;
    T2[i][0]  = - Psi[i][1];
    T2[i][1]  = 0;
    T3[i][0]  = Psi[i][0];
    T3[i][1]  = 0;
  }

  /*-- compute:   T1 = v_x,  T2 = v_y,  T3 = laplacian(u) --*/
  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);
  fft_laplacian(T3, grid);

  /*--  compute dv = laplacian(u) - b ( grad(v) )^2 / 2 --*/ 
  /*--  compute [T1,T2] = (1 + a u) grad(v) --*/

  for (i=0; i<N; i++){

    Psi[i][1] = - T3[i][0] + 0.5 * b * (T1[i][0]*T1[i][0] + T2[i][0]*T2[i][0]);

    T1[i][0]  = (1 + a * Psi[i][0]) * T1[i][0];
    T1[i][1]  = 0;
    
    T2[i][0]  = (1 + a * Psi[i][0]) * T2[i][0];
    T2[i][1]  = 0;
  }

  /*--  compute du = - div( (1 + a u) grad(v) )  --*/

  
  fft_deriv_x(T1, grid);
  fft_deriv_y(T2, grid);

  for (i=0; i<N; i++){

    Psi[i][0] = - T1[i][0] - T2[i][0];

  }

  
}

/* ---------------------------------------------------------------- */
