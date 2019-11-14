#include "header.h"

static  fftw_complex  *S, *T1, *T2;

static  double         coefRho, coefNu;

static  int N0;
static  int np;

/* ---------------------------------------------------------------- */

void rhs_init(geom_ptr geom, phys_ptr phys, int np_in)
{

  int Ngrids, N;

  np = np_in;

  Ngrids =  geom->Ngrids;
  N0 =  geom->N0;

  N = N0*pow(2, Ngrids-1);
  N = N*N/np;
    
  T1 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));
  T2 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));

  coefRho = phys->coefRho;
  coefNu  = phys->coefNu;

}

/* ---------------------------------------------------------------- */


void rhs_init_S(fftw_complex *s)
{

  S = s;

}

/* ---------------------------------------------------------------- */


void rhs_compute(int grid)   // S = rhs(S) 
{

  int i, N;
  double u, v, pp;

  N = N0 * pow(2,grid);
  N = N*N/np;

  for (i=0; i<N; i++){
    u = S[i][0];
    v = S[i][1];
    T1[i][0]  = u;
    T1[i][1]  = v;
    T2[i][0]  = u*u + v*v;
    T2[i][1]  = 0;
  }

  fft_laplacian(T1, grid);
  fft_davey_stewartson_phi_x(T2, grid, coefNu);

  for (i=0; i<N; i++){

    u    = S[i][0];
    v    = S[i][1];
    pp   = u*u + v*v;

    S[i][0] = - T1[i][1] - pp*v  + coefRho * (u*T2[i][1] + v*T2[i][0]);
    S[i][1] =   T1[i][0] + pp*u  - coefRho * (u*T2[i][0] - v*T2[i][1]);

  }

}

/* ---------------------------------------------------------------- */
