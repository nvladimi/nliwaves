#include "header.h"


static  fftw_complex  *Psi, *S, *S0;

extern void   rhs_compute(int grid);
extern void   rhs_init_S(fftw_complex *s);

static void   rk_zero_out(int N);
static void   rk_add_slope(int N, int w);
static void   rk_part_step(int N, double dt);
static void   rk_full_step(int N, double dt);

static int N0;
static int np;

/*--------------------------------------------------------------*/

void rk_init(geom_ptr geom, fftw_complex *psi, int np_in)
{
  int N, Ngrids;

  np = np_in;

  Ngrids =  geom->Ngrids;
  N0 =  geom->N0;

  N = N0*pow(2, Ngrids-1);
  N = N*N/np;

  S  =  (fftw_complex *) malloc( N * sizeof(fftw_complex));
  S0 =  (fftw_complex *) malloc( N * sizeof(fftw_complex));

  Psi = psi;

  rhs_init_S(S);
}

/*--------------------------------------------------------------*/

void rk_one_step(int grid, double dt)   // Uo := Uo + dU
{

  int N;
  
  N = N0 * pow(2,grid);
  N = N*N/np;

  rk_zero_out(N);

  /*-- stage 1 --*/

  rk_part_step(N, 0);
  rhs_compute(grid);
  rk_add_slope(N, 1);


  /*-- stage 2 --*/

  rk_part_step(N, dt/2);
  rhs_compute(grid);
  rk_add_slope(N, 2);


  /*-- stage 3 --*/

  rk_part_step(N, dt/2);
  rhs_compute(grid);
  rk_add_slope(N, 2);


  /*-- stage 4 --*/

  rk_part_step(N, dt);
  rhs_compute(grid);
  rk_add_slope(N, 1);


 /*-- update solution --*/

  rk_full_step(N, dt/6);

  
}

/*--------------------------------------------------------------*/

void rk_zero_out(int N)              // S=0,  S0=0;
{
  memset(S,  0, N*sizeof(fftw_complex));
  memset(S0, 0, N*sizeof(fftw_complex));
}

void rk_part_step(int N, double dt)  // U = Uo + S*h,   S: slope -> solution
{
  int i;
  for (i=0; i<N; i++) { 
    S[i][0] = Psi[i][0] + S[i][0]*dt; 
    S[i][1] = Psi[i][1] + S[i][1]*dt; 
  }
}

//void  compute_rhs(int N)          // S: solution -> slope 

void rk_add_slope(int N, int w)     // S0 = S0 + S*w
{
  int i;
  for (i=0; i<N; i++) { 
    S0[i][0] += S[i][0] * w; 
    S0[i][1] += S[i][1] * w; 
  }
}

void rk_full_step(int N, double dt)  // Uo = Uo + So*h, 
{
  int i;
  for (i=0; i<N; i++) { 
    Psi[i][0] += S0[i][0]*dt; 
    Psi[i][1] += S0[i][1]*dt; 
  }
}

/*--------------------------------------------------------------*/
